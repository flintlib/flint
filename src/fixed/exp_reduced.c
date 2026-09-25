/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "fixed.h"

/* fixed_exp_reduced: exp(t) of a reduced argument t < 2^-r,
   r >= 16 (the series algorithms 1 and 2 require r >= 32), into y
   (wn + 1 limbs: wn fraction limbs and a units limb), independent
   of any particular argument-reduction scheme.

   alg selects the method:
       0  the tuned automatic choice (see the selection logic at the
          bottom and tune/tune-exp-reduced.c)
       1  the direct exponential rectangular-splitting series
       2  the sinh rectangular-splitting series plus a squaring and
          a square root
       3  one bit-burst step -- t = x1 + x2 with x1 the leading
          slice down to the doubled limb boundary and x2 below it --
          combining exp(x1) from the mpn binary splitting with
          exp(x2) from this function again at the doubled rate
       4  the full bit-burst algorithm: exp(t) = prod_k exp(x_k)
          with slices doubling in length on limb boundaries, each
          factor from the binary splitting; asymptotically
          quasi-optimal for very large wn (or very large r)

   The slice boundaries sit on limbs: the leading slice runs from
   below the r zero bits down to limb wn - D, D = 2 max(r/64, 1),
   so its value is below 2^-r and the remainder below 2^-64D; the
   binary splitting reads the slice as an integer scaled by
   2^(-64 D), so an r not divisible by FLINT_BITS merely leaves a
   few zero bits at the top of the slice.  All the arithmetic
   around the series kernels and the exact binary splitting is
   fball ball arithmetic (see _fball_exp_reduced), whose rigorous
   bound is checked against FIXED_EXP_REDUCED_MAX_ERR on export (it
   comes to about 27 ulps for the one-step cascade, 1 for the full
   burst).  An earlier driver in hand-written mpn arithmetic
   (windowed middle products, limb-aligned frames and an informal
   error accounting) measured the same speed at every size and was
   dropped. */

/* Use the sinh series (half the terms) plus a squaring and a square
   root once the direct exp series gets long enough.  Measured
   crossovers on x86-64: for r < 64 (series through the pre32 path)
   sinh wins from about 45 terms (n ~ 24 at r = 32); for r >= 64 the
   windowed pre64 series is more efficient per term and flint_mpn_sqrtrem
   sets a higher floor, moving the crossover to about 128 terms
   (n ~ 2r).  The margins near the crossovers are within a few
   percent, so the exact placement is not critical. */
/* sinh from this many series terms (64 wn / r): the crossover rises
   with r (tune-exp-reduced on the target machine: (64, 128] up to
   r = 256, (171, 256] for r = 384 .. 2048, (256, 512] at 4096) */
#define EXP_USE_SINH(wn, r) \
    (FLINT_BITS * (wn) >= (((r) >= 3072) ? 300 : ((r) >= 384) ? 200 \
        : ((r) >= 64) ? 128 : 45) * (slong) (r))

/* When the direct series would need many terms, evaluate sinh(t)
   instead -- the odd series has half the terms -- and reconstruct
   exp(t) = sinh(t) + sqrt(1 + sinh(t)^2), a ball squaring, sum and
   square root on the imported sinh (FIXED_SINH_RS_MAX_ERR = 15 ulp),
   which comes to about 17 ulps.  One bit-burst step doubles the
   convergence rate: for t < 2^-r with r a whole number of limbs,
   split t = x1 + x2 where x1 is the leading r-bit chunk (bits
   [r, 2r)) and x2 < 2^-2r the remainder, and compute
   exp(t) = exp(x1) exp(x2), exp(x1) = (Q B^QE + T) / (Q B^QE) from
   the mpn binary splitting _fixed_exp_sum_bs_powtab -- the chunk is
   an r-bit rational, so the splitting tree stays small -- and
   exp(x2) from this function at the doubled rate (possibly bursting
   again). */
/* terms threshold for one burst step; the slice mechanics work for
   any r >= 32 since the boundaries sit on limbs.  Measured
   crossovers (tune-exp-reduced on the target machine): (256, 512]
   at r = 64 and 128, (341, 683] .. (512, 1024] for r = 192 .. 1024,
   (683, 1365] from r = 1536. */
#define FIXED_EXP_BURST_TERMS(r) \
    (((r) >= 1536) ? 1000 : ((r) > 128) ? 600 : 512)
#define EXP_USE_BURST(wn, r) \
    (FLINT_BITS * (wn) >= FIXED_EXP_BURST_TERMS(r) * (slong) (r))

/* Full bit-burst from this many terms.  The development VM's value
   (4096) sat below the crossover measured on the target machine at
   every r: (5461, 8192] for r >= 64, (8192, 16384] at r = 32,
   (16384, 32768] at r = 16 -- the one-step cascade, which finishes
   by rectangular splitting, beats the full burst for longer than
   the earlier tuning found.  (The old rationale: one burst step
   from ~512 terms, later than under the per-level divisions whose
   small level-0 divisor was cheaper than the generic rational
   accumulation for one slice.) */
#define FIXED_EXP_FULLBURST_TERMS(r) \
    (((r) >= 64) ? 8192 : ((r) >= 32) ? 16384 : 32768)

/* Bit length of an unsigned mpn (x, xn), xn >= 1: the position of
   its most significant set bit, i.e. 64 (xn - 1) + bitcount(top).
   For a normalized (top limb nonzero) operand this is its exact
   floor(log2) + 1. */
static slong
mpn_bits(nn_srcptr x, slong xn)
{
    return FLINT_BITS * (xn - 1) + FLINT_BIT_COUNT(x[xn - 1]);
}

#ifndef EXP_BURST_GUARD
#define EXP_BURST_GUARD 3
#endif

/* ==== the algorithms in fball arithmetic ================================

   exp(t) as a ball, with the exponent and error bookkeeping left to
   fball: only the two series kernels are imported (fixed_exp_rs and
   fixed_sinh_rs, with their documented one-sided errors), the sinh
   reconstruction sinh + sqrt(1 + sinh^2) is a ball squaring, sum and
   square root, and in the burst each slice's exact splitting output
   T / (Q B^QEk) becomes the factor Q B^QEk + T by one ball addition
   (which reads only the top cap + 2 limbs of either operand), the
   products NUM *= factor and DEN *= Q are ball products at cap limbs
   (their operands truncated to the precision that contributes), the
   series remainder of the one-step variant is this function again
   at the doubled rate (the cascade), and the finish is one ball
   division.  The bound of the result is rigorous, and checked
   against the documented budget on export (it lands far below). */
static void
_fball_exp_reduced(fball_t res, nn_srcptr t, slong wn, flint_bitcnt_t r,
    int alg)
{
    slong cap = wn + EXP_BURST_GUARD;

    if (alg == 0)
    {
        if (FLINT_BITS * (ulong) wn >= (ulong) FIXED_EXP_FULLBURST_TERMS(r) * r)
            alg = 4;
        else if (r < 32 || EXP_USE_BURST(wn, r))
            alg = 3;
        else if (EXP_USE_SINH(wn, r))
            alg = 2;
        else
            alg = 1;
    }

    if (alg == 1)
    {
        nn_ptr y;
        TMP_INIT;
        TMP_START;
        y = TMP_ALLOC((wn + 1) * sizeof(ulong));
        fixed_exp_rs(y, t, wn);
        fball_set_mpn_2exp(res, y, wn + 1, -FLINT_BITS * wn);
        fball_add_error(res, FIXED_EXP_RS_MAX_ERR(wn), -wn);
        TMP_END;
    }
    else if (alg == 2)
    {
        /* exp = sinh + sqrt(1 + sinh^2) */
        fball_t s, u;
        nn_ptr y;
        TMP_INIT;
        TMP_START;
        y = TMP_ALLOC((wn + 1) * sizeof(ulong));
        fball_init(s);
        fball_init(u);
        fixed_sinh_rs(y, t, wn);
        fball_set_mpn_2exp(s, y, wn + 1, -FLINT_BITS * wn);
        fball_add_error(s, FIXED_SINH_RS_MAX_ERR(wn), -wn);
        fball_mul(u, s, s, cap);
        fball_set_ui(res, 1);
        fball_add(u, u, res, cap);
        fball_sqrt(u, u, cap);
        fball_add(res, u, s, cap);
        fball_clear(s);
        fball_clear(u);
        TMP_END;
    }
    else
    {
        int levels = (alg == 4) ? 0 : 1;
        slong L[FLINT_BITS + 2];
        slong nb = 0, k;
        fball_t num, den, f, fq;

        L[nb++] = FLINT_MAX((slong) r / FLINT_BITS, 1);
        while (L[nb - 1] < wn)
        {
            slong nxt = FLINT_MIN(2 * L[nb - 1], wn);
            if (levels > 0 && nb > levels)
                nxt = wn;
            L[nb] = nxt;
            nb++;
        }
        nb--;
        if (nb == 0)
        {
            L[1] = wn;
            nb = 1;
        }

        fball_init(num);
        fball_init(den);
        fball_init(f);
        fball_init(fq);
        fball_set_ui(num, 1);
        fball_set_ui(den, 1);

        for (k = nb - 1; k >= 0; k--)
        {
            int series = (levels > 0 && k >= levels);
            TMP_INIT;

            TMP_START;

            if (series)
            {
                /* the residual below limb depth L[k], at the doubled
                   rate by the tuned choice */
                nn_ptr xres;

                xres = TMP_ALLOC(wn * sizeof(ulong));
                flint_mpn_copyi(xres, t, wn);
                flint_mpn_zero(xres + wn - L[k], L[k]);
                _fball_exp_reduced(f, xres, wn,
                    (flint_bitcnt_t) (FLINT_BITS * L[k]), 0);
            }
            else
            {
                slong D = L[k + 1];
                nn_srcptr u = t + (wn - D);
                slong xn = D - (k ? L[k] : 0);
                slong N, tn, qn, QEk;
                nn_ptr T, Q;

                while (xn > 0 && u[xn - 1] == 0)
                    xn--;
                if (xn == 0)
                {
                    TMP_END;
                    continue;
                }
                while (xn > 1 && u[0] == 0)
                {
                    u++;
                    xn--;
                    D--;
                }

                N = _fixed_exp_bs_num_terms((flint_bitcnt_t)
                    (FLINT_BITS * D - mpn_bits(u, xn)),
                    FLINT_BITS * wn + 64);

                {
                    slong t_alloc = N * (D + 2) + 4;
                    slong q_alloc2 = (N * FLINT_BIT_COUNT((ulong) N + 1))
                        / FLINT_BITS + 3;

                    T = TMP_ALLOC((t_alloc + q_alloc2) * sizeof(ulong));
                    Q = T + t_alloc;
                }

                _fixed_exp_sum_bs_powtab(T, &tn, Q, &qn, &QEk,
                    u, xn, D, N);

                /* the factor Q B^QEk + T over the denominator Q B^QEk */
                fball_set_mpn_2exp(fq, Q, qn, FLINT_BITS * QEk);
                fball_set_mpn_2exp(f, T, tn, 0);
                fball_add(f, f, fq, cap);
                fball_mul(den, den, fq, cap);
            }

            fball_mul(num, num, f, cap);
            TMP_END;
        }

        fball_div(res, num, den, wn + 2);

        fball_clear(num);
        fball_clear(den);
        fball_clear(f);
        fball_clear(fq);
    }
}

/* the fball algorithms exported in the fixed-point format */
static void
_fixed_exp_reduced_fball(nn_ptr y, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg)
{
    fball_t v, one;
    double bound;

    fball_init(v);
    fball_init(one);
    _fball_exp_reduced(v, t, wn, r, alg);
    /* exp(t) - 1 in [0, 1) */
    fball_set_ui(one, 1);
    fball_sub(v, v, one, wn + 2);
    bound = fball_get_fixed(y, wn, v);
    if (!(bound < FIXED_EXP_REDUCED_MAX_ERR))
        flint_throw(FLINT_ERROR, "fixed_exp_reduced (fball): error bound %g ulps\n", bound);
    y[wn] = 1;
    fball_clear(v);
    fball_clear(one);
}

void
fixed_exp_reduced(nn_ptr y, nn_srcptr t, slong wn, flint_bitcnt_t r,
    int alg)
{
    FLINT_ASSERT(r >= 16);
    FLINT_ASSERT(alg == 0 || alg >= 3 || r >= 32);
    _fixed_exp_reduced_fball(y, t, wn, r, alg);
}
