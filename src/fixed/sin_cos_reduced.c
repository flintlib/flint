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

/* fixed_sin_cos_reduced: sin(t) and g = 1 - cos(t) of a reduced
   argument t < 2^-r, r >= 16 (the series algorithms 1 and 2 require
   r >= 32), into (ysin, wn) and (yg, wn) -- wn fraction limbs each;
   sin(t) < 2^-r and g < 2^(-2r-1), so neither carries a units limb,
   though both buffers must have room for wn + 1 limbs (the top limb
   is scratch).  This is the trigonometric analogue of
   fixed_exp_reduced, independent of any particular
   argument-reduction scheme (the tangent half-angle reconstruction
   consumes exactly this pair).

   alg selects the method:
       0  the tuned automatic choice
       1  the direct sine + cosine rectangular-splitting series
       2  the sine rectangular-splitting series plus a squaring and
          a square root (the "sin sqrt" trick: the odd series has
          half the terms, and g = 1 - sqrt(1 - sin^2))
       3  one bit-burst step -- t = x1 + x2 with x1 the leading
          slice down to the tripled first limb boundary and x2 below
          it -- combining exp(i x1) from the Gaussian mpn binary
          splitting with exp(i x2) from the tuned series at the
          raised rate
       4  the full bit-burst algorithm: exp(i t) = prod_k exp(i x_k)
          with slice lengths TRIPLING on limb boundaries, each
          factor from the Gaussian binary splitting; asymptotically
          quasi-optimal for very large wn (or very large r)

   The burst mirrors fixed_exp_reduced at LIMB granularity: the
   boundary ladder is a list of limb counts (tripling), a slice is a
   pointer + length into t with a whole-limb mask, dead low limbs of
   a slice fold into its frame (u B^-D = (u/B^z) B^-(D-z)), and each
   slice's factor is the exact Gaussian rational

       exp(i x_k) ~ (Q_k B^{QE_k} - A_k + i B_k) / (Q_k B^{QE_k})

   truncated at N_k series terms, all exponents limb counts.
   Instead of combining per level, the loop accumulates the complex
   numerator NUM = prod (Q_k B^{QE_k} - A_k + i B_k) and the real
   denominator DEN = prod Q_k as balls at wn + 3 limbs (one
   transform-sharing complex product per slice, fball_mul_complex),
   deepest slice first so every accumulation product is balanced
   against the content gathered so far, and finishes with TWO ball
   divisions (sine and cosine) against the single accumulated
   denominator -- where arb_sin_cos_arf_bb spends one full-precision
   square root PER SLICE (cos_k from sin_k) plus a full-precision
   complex product tree.  g comes from the cosine by subtraction:
   the ~2r-bit cancellation is harmless because g is only needed to
   absolute 2^(-64 wn), like the direct series path.  The bounds are
   rigorous, checked against FIXED_SIN_COS_REDUCED_MAX_ERR = 96 ulps
   on export (they come to a few dozen).  An earlier driver in
   hand-written mpn arithmetic (windowed middle products, limb
   frames, flint_mpn_mulhigh_n_complex, an informal error accounting)
   measured the same speed phase for phase and was dropped. */

#define TRIG_BURST_GUARD 3

/* Thresholds measured on the development VM (terms = 64 wn / r;
   sweeps in dev/notes), after the per-slice square-root hybrid and
   the truncated trees.  The sqrt-trick crossover carries over the
   50-term value the tangent half-angle code was tuned to (measured
   ~24-48 here, margins of a few percent).  The one-burst-step
   crossover is strongly r-dependent: at r < 64 a single step
   halves a LONG series and wins from ~256-384 terms, while for
   r >= 64 the splitting integers' fat terms push it to ~2048-3072.
   The one-step variant CASCADES: its series remainder runs at the
   tripled rate and redispatches, so below the full-burst threshold
   the effective algorithm is a few burst slices finished by
   rectangular splitting once the remaining series is short -- and
   that hybrid beats the full bit-burst (and arb_sin_cos_arf_bb) by
   1.2-1.3x in the 1e5-4e6 bit range.  The full burst only draws
   level from ~4e5 series terms (~1e7 bits at r = 24; measured at
   r = 24 and r = 32, the terms crossover agreeing), so a single
   high threshold covers all r; the tuning is self-referential
   (each threshold change alters the cascade being measured) and
   these values are the fixed point of that iteration on the
   development VM.

   Re-swept for r >= 64 on the fft_small build over n = 256 .. 15625
   limbs and r = 32 .. 768 (the diophantine (multi-prime) reduction of
   sin_cos_diophantine.c lands at r ~ 100-300): the one-step crossover
   in terms/r grows slowly with n -- about 500 at 512 limbs, 800 at
   1024, 1300 at 2677, 1500 at 6000, below 1950 at 15625 -- so the
   earlier 2560 sat above it at every size from 1024 limbs up,
   choosing the sine+sqrt path in a band where the burst step is
   5-25% faster.  1500 is the single constant that makes the right
   call at 2677, 6000 and 15625 limbs and costs about 6% in one cell
   at 1024 limbs, r = 96; it brings sin_cos_diophantine's default at
   640k bits from 57.8 to 45.2 ms and its 13-prime configuration at
   171k bits from 13.0 to 10.9 ms.  The target machine's
   tune-sin-cos-reduced (wn = 16 .. 131072, r = 16 .. 4096) puts the
   same crossover in (1365, 2048] for every r >= 64, consistent with
   1500, and (256, 512] at r = 32, consistent with 320. */
/* sine + sqrt from this many terms: 50 up to r = 1536; the
   crossover moves to (64, 128] at r = 2048 .. 3072 and (128, 256]
   at r = 4096 (tune-sin-cos-reduced on the target machine) */
#define FIXED_TRIG_REDUCED_SINSQRT_TERMS(r) \
    (((r) >= 4096) ? 200 : ((r) >= 2048) ? 100 : 50)
#ifndef FIXED_TRIG_BURST_TERMS_SMALL_R
#define FIXED_TRIG_BURST_TERMS_SMALL_R 320
#endif
#ifndef FIXED_TRIG_BURST_TERMS
#define FIXED_TRIG_BURST_TERMS 1500
#endif
/* the full bit-burst never won on the target machine up to 524288
   terms per r (8.4 10^6 bits at r = 16); the one-step cascade covers
   that range, so the switch sits above it */
#ifndef FIXED_TRIG_FULLBURST_TERMS
#define FIXED_TRIG_FULLBURST_TERMS 1048576
#endif
/* Per-slice choice inside the burst: from this many series terms
   the slice's 1 - cos track (two of the four heavy tree
   multiplications per merge) costs more than recovering the cosine
   window by a fixed-point square root -- one unbalanced division
   by the slice denominator, a squaring, the square root and one
   short product against Q.  Below it the joint tree wins (the
   square-root path's costs are flat in N). */
#ifndef FIXED_TRIG_SLICE_SQRT_TERMS
#define FIXED_TRIG_SLICE_SQRT_TERMS 64
#endif
#define TRIG_USE_SINSQRT(wn, r) \
    (FLINT_BITS * (wn) >= FIXED_TRIG_REDUCED_SINSQRT_TERMS(r) * (slong) (r))
#define TRIG_USE_BURST(wn, r) \
    (FLINT_BITS * (wn) >= (((r) < 64) ? FIXED_TRIG_BURST_TERMS_SMALL_R \
        : FIXED_TRIG_BURST_TERMS) * (slong) (r))


/* ==== the algorithms in fball arithmetic ================================

   sin(t) and cos(t) as balls (the export forms g = 1 - cos): only the
   series kernels are imported with their documented errors, the sine
   + square root reconstruction is a ball squaring, sum and square
   root, and in the burst each slice's exact Gaussian splitting output
   becomes the factor FC + i FS by ball additions (FC = Q B^QE - A, or
   sqrt((Q B^QE)^2 - FS^2) for the sine-only slices) and one short
   product (FS = x (Q B^QE - B) B^-D), the complex accumulation
   NUM *= FC + i FS is four ball products and two sums at cap limbs,
   DEN *= Q one product, the series remainder of the one-step variant
   is this function again at the tripled rate, and the finish is two
   ball divisions.  The bounds are rigorous and checked against the
   documented budget on export. */
static void
_fball_sin_cos_reduced(fball_t rs, fball_t rc, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg)
{
    slong cap = wn + TRIG_BURST_GUARD;

    if (alg == 0)
    {
        ulong terms = FLINT_BITS * (ulong) wn;
        if (terms >= (ulong) FIXED_TRIG_FULLBURST_TERMS * r)
            alg = 4;
        else if (r < 32 || TRIG_USE_BURST(wn, r))
            alg = 3;
        else if (TRIG_USE_SINSQRT(wn, r))
            alg = 2;
        else
            alg = 1;
    }

    if (alg == 1)
    {
        nn_ptr ss, cc;
        TMP_INIT;
        TMP_START;
        ss = TMP_ALLOC(2 * (wn + 2) * sizeof(ulong));
        cc = ss + (wn + 2);
        fixed_sin_cos_rs(ss, cc, t, wn);
        fball_set_mpn_2exp(rs, ss, wn + 1, -FLINT_BITS * wn);
        fball_add_error(rs, FIXED_SIN_COS_RS_MAX_ERR(wn), -wn);
        fball_set_mpn_2exp(rc, cc, wn + 1, -FLINT_BITS * wn);
        fball_add_error(rc, FIXED_SIN_COS_RS_MAX_ERR(wn), -wn);
        TMP_END;
    }
    else if (alg == 2)
    {
        /* cos = sqrt(1 - sin^2) */
        fball_t u;
        nn_ptr ss;
        TMP_INIT;
        TMP_START;
        ss = TMP_ALLOC((wn + 2) * sizeof(ulong));
        fball_init(u);
        fixed_sin_rs(ss, t, wn);
        fball_set_mpn_2exp(rs, ss, wn + 1, -FLINT_BITS * wn);
        fball_add_error(rs, FIXED_SIN_RS_MAX_ERR(wn), -wn);
        fball_mul(u, rs, rs, cap);
        fball_set_ui(rc, 1);
        fball_sub(u, rc, u, cap);
        fball_sqrt(rc, u, cap);
        fball_clear(u);
        TMP_END;
    }
    else
    {
        int levels = (alg == 4) ? 0 : 1;
        slong L[FLINT_BITS + 2];
        slong nb = 0, k;
        fball_t nc, ns, den, fc, fs, fq, u, v;

        L[nb++] = FLINT_MAX((slong) r / FLINT_BITS, 1);
        while (L[nb - 1] < wn)
        {
            slong nxt = FLINT_MIN(3 * L[nb - 1], wn);
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

        fball_init(nc);
        fball_init(ns);
        fball_init(den);
        fball_init(fc);
        fball_init(fs);
        fball_init(fq);
        fball_init(u);
        fball_init(v);
        fball_set_ui(nc, 1);
        fball_set_ui(den, 1);

        for (k = nb - 1; k >= 0; k--)
        {
            int series = (levels > 0 && k >= levels);
            TMP_INIT;

            TMP_START;

            if (series)
            {
                nn_ptr xres;

                xres = TMP_ALLOC(wn * sizeof(ulong));
                flint_mpn_copyi(xres, t, wn);
                flint_mpn_zero(xres + wn - L[k], L[k]);
                _fball_sin_cos_reduced(fs, fc, xres, wn,
                    (flint_bitcnt_t) (FLINT_BITS * L[k]), 0);
            }
            else
            {
                slong D = L[k + 1];
                nn_srcptr x = t + (wn - D);
                slong xn = D - (k ? L[k] : 0);
                slong N, an, bn2, qn, ae, be, QEk;
                nn_ptr A, B, Q;
                int slice_sqrt;

                while (xn > 0 && x[xn - 1] == 0)
                    xn--;
                if (xn == 0)
                {
                    TMP_END;
                    continue;
                }
                while (xn > 1 && x[0] == 0)
                {
                    x++;
                    xn--;
                    D--;
                }

                {
                    slong ubits = FLINT_BITS * (xn - 1)
                        + FLINT_BIT_COUNT(x[xn - 1]);
                    N = _fixed_exp_bs_num_terms(
                        (flint_bitcnt_t) (FLINT_BITS * D - ubits),
                        FLINT_BITS * wn + 64);
                    N = FLINT_MAX(1, N / 2);
                    if (N > 10000)
                        while (N % 128 != 0)
                            N++;
                    if (N > 1000)
                        while (N % 16 != 0)
                            N++;
                    if (N > 100)
                        while (N % 2 != 0)
                            N++;
                }

                {
                    slong qb2 = (N * 2
                        * FLINT_BIT_COUNT(2 * (ulong) N + 1))
                        / FLINT_BITS + 3;
                    A = TMP_ALLOC((2 * (cap + 2 + qb2 + 8) + qb2)
                        * sizeof(ulong));
                    B = A + (cap + 2 + qb2 + 8);
                    Q = B + (cap + 2 + qb2 + 8);
                }

                slice_sqrt = (N >= FIXED_TRIG_SLICE_SQRT_TERMS);
                _fixed_sin_cos_sum_bs_powtab(
                    slice_sqrt ? NULL : A, &an, &ae, B, &bn2, &be,
                    Q, &qn, &QEk, x, xn, D, N, cap + 2);

                /* the denominator Q B^QEk of the slice */
                fball_set_mpn_2exp(fq, Q, qn, FLINT_BITS * QEk);
                fball_mul(den, den, fq, cap);

                /* FS = x (Q B^QEk - B B^be) B^-D */
                fball_set_mpn_2exp(fs, B, bn2, FLINT_BITS * be);
                fball_sub(fs, fq, fs, cap);
                fball_set_mpn_2exp(u, x, xn, -FLINT_BITS * D);
                fball_mul(fs, fs, u, cap);

                if (!slice_sqrt)
                {
                    /* FC = Q B^QEk - A B^ae */
                    fball_set_mpn_2exp(fc, A, an, FLINT_BITS * ae);
                    fball_sub(fc, fq, fc, cap);
                }
                else
                {
                    /* FC = sqrt((Q B^QEk)^2 - FS^2) */
                    fball_mul(u, fq, fq, cap);
                    fball_mul(v, fs, fs, cap);
                    fball_sub(u, u, v, cap);
                    fball_sqrt(fc, u, cap);
                }
            }

            /* NUM *= FC + i FS */
            if (nc->size == 1 && nc->d[0] == 1 && nc->err == 0
                && nc->exp == 1 && ns->size == 0 && ns->err == 0)
            {
                fball_set(nc, fc);
                fball_set(ns, fs);
            }
            else
                fball_mul_complex(nc, ns, nc, ns, fc, fs, cap);
            TMP_END;
        }

        fball_div(rs, ns, den, wn + 2);
        fball_div(rc, nc, den, wn + 2);

        fball_clear(nc);
        fball_clear(ns);
        fball_clear(den);
        fball_clear(fc);
        fball_clear(fs);
        fball_clear(fq);
        fball_clear(u);
        fball_clear(v);
    }
}

/* the fball algorithms exported in the fixed-point format */
static void
_fixed_sin_cos_reduced_fball(nn_ptr ysin, nn_ptr yg, nn_srcptr t,
    slong wn, flint_bitcnt_t r, int alg)
{
    fball_t s, c, one;
    double bound;

    fball_init(s);
    fball_init(c);
    fball_init(one);
    _fball_sin_cos_reduced(s, c, t, wn, r, alg);
    bound = fball_get_fixed(ysin, wn, s);
    if (!(bound < FIXED_SIN_COS_REDUCED_MAX_ERR))
        flint_throw(FLINT_ERROR, "fixed_sin_cos_reduced (fball): sine error bound %g ulps\n", bound);
    /* g = 1 - cos in [0, 1) */
    fball_set_ui(one, 1);
    fball_sub(c, one, c, wn + 2);
    bound = fball_get_fixed(yg, wn, c);
    if (!(bound < FIXED_SIN_COS_REDUCED_MAX_ERR))
        flint_throw(FLINT_ERROR, "fixed_sin_cos_reduced (fball): cosine error bound %g ulps\n", bound);
    fball_clear(s);
    fball_clear(c);
    fball_clear(one);
}

void
fixed_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg)
{
    FLINT_ASSERT(r >= 16);
    FLINT_ASSERT(alg == 0 || alg >= 3 || r >= 32);
    _fixed_sin_cos_reduced_fball(ysin, yg, t, wn, r, alg);
}
