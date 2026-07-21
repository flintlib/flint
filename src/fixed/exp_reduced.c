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
          exp(x2) from the sinh series at the doubled rate
       4  the full bit-burst algorithm: exp(t) = prod_k exp(x_k)
          with slices doubling in length on limb boundaries, each
          factor from the binary splitting; asymptotically
          quasi-optimal for very large wn (or very large r)

   The slice boundaries sit on limbs: the leading slice runs from
   below the r zero bits down to limb wn - D, D = 2 max(r/64, 1),
   so its value is below 2^-r and the remainder below 2^-64D; the
   binary splitting reads the slice as an integer scaled by
   2^(-64 D), so an r not divisible by FLINT_BITS merely leaves a
   few zero bits at the top of the slice.  Error: the series paths
   carry their documented budgets (~20 ulp for sinh); each burst
   level adds at most ~4 ulps (a floored short division within 2 and
   a middle product taken two limbs below the kept window within 2),
   and alg 4 uses at most log2(wn) levels. */

/* Use the sinh series (half the terms) plus a squaring and a square
   root once the direct exp series gets long enough.  Measured
   crossovers on x86-64: for r < 64 (series through the pre32 path)
   sinh wins from about 45 terms (n ~ 24 at r = 32); for r >= 64 the
   windowed pre64 series is more efficient per term and mpn_sqrtrem
   sets a higher floor, moving the crossover to about 128 terms
   (n ~ 2r).  The margins near the crossovers are within a few
   percent, so the exact placement is not critical. */
/* working precision from which the sinh-path square root runs
   through fixed_sqrt_newton instead of mpn_sqrtrem; kept above the
   measured range on the development machine (see the comment at the
   call site), lower it to experiment */
#ifndef FIXED_EXP_SQRT_NEWTON_CUTOFF
#define FIXED_EXP_SQRT_NEWTON_CUTOFF 2000
#endif

#define EXP_USE_SINH(wn, r) \
    (FLINT_BITS * (wn) >= (((r) >= 64) ? 128 : 45) * (slong) (r))

/* exp(t) of the reduced argument t < 2^-r into y (wn + 1 limbs:
   wn fraction limbs and a units limb); the series functions pick the
   32- or 64-bit internal range from the top limb of t.  When the
   direct series would need many terms, evaluate sinh(t) instead --
   the odd series has half the terms -- and reconstruct
   exp(t) = sinh(t) + sqrt(1 + sinh(t)^2), costing one squaring and
   one square root.  The sinh error (FIXED_SINH_RS_MAX_ERR = 15 ulp),
   the truncated squaring and the floored integer square root
   together contribute some 20 ulp, amplified below e by the
   reconstruction: this is part of the constant term of
   FIXED_EXP_BITWISE_RS_MAX_ERR. */
/* One bit-burst step doubles the convergence rate: for t < 2^-r
   with r a whole number of limbs, split t = x1 + x2 where x1 is the
   leading r-bit chunk (bits [r, 2r)) and x2 < 2^-2r the remainder,
   and compute exp(t) = exp(x1) exp(x2).  exp(x1) = 1 + s1 comes
   from the mpn binary splitting _fixed_exp_sum_bs_powtab -- the
   chunk is an r-bit rational, so the splitting tree stays small --
   and exp(x2) = 1 + f2 recurses into this function at the doubled
   rate (possibly bursting again).  The combination
   (1 + s1)(1 + f2) = 1 + s1 + f2 + s1 f2 costs two additions and
   one balanced middle product for s1 f2 < 2^-3r, taken two limbs
   below the kept window per the mulmid deficit bound.  s1 is a
   floored quotient within 2 ulps (shift truncation plus tdiv), the
   middle product adds at most 2 more; all far inside the sinh
   reconstruction budget. */
#ifndef FIXED_EXP_BURST_TERMS
#define FIXED_EXP_BURST_TERMS 512
#endif
/* terms threshold for one burst step; the slice mechanics work for
   any r >= 32 since the boundaries sit on limbs */
#define EXP_USE_BURST(wn, r) \
    (FLINT_BITS * (wn) >= FIXED_EXP_BURST_TERMS * (slong) (r))

/* Thresholds from tune-exp-reduced on the development VM after the
   single-division restructure (terms = 64 wn / r): one burst step
   from ~512 terms -- LATER than under the per-level divisions,
   whose small level-0 divisor was cheaper than the generic rational
   accumulation for one slice -- and the full bit-burst from ~4096
   terms, EARLIER than before since the cascade is what the
   restructure optimizes.  Retune on target hardware. */
#ifndef FIXED_EXP_FULLBURST_TERMS
#define FIXED_EXP_FULLBURST_TERMS 4096
#endif

/* Assemble the bit-burst factor F = Q 2^Qexp + T as
   value = (F, *fn) 2^(*fexp): exactly when it fits within
   cap + 2 limbs (then *fexp = 0), and otherwise as its top window
   of about cap + 1 limbs -- Q shifted into place and T's
   limb-aligned tail above the window bottom added -- with the
   dropped bits, under one bit at the limb-aligned bottom, going
   into *fexp.  F must have room for cap + 3 limbs.  Exposed
   non-statically so that the window arithmetic, which the burst
   only reaches in specific precision regimes, can be unit-tested
   directly (t-exp_sum_bs). */
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

/* ==== the bit-burst driver ==============================================

   Slice boundaries snap to LIMBS and every
   power-of-two frame (fexp, nexp, dexp, QE, the final placement) is
   a LIMB count.  Consequences, all following from limb granularity:

   - slice extraction is a pointer and a length: bits shallower than
     a boundary live in separate limbs, so there is no rshift and no
     masking (the top slice's leading bits above depth r are zero in
     t itself);
   - the series-leaf residual clears whole limbs, no partial-limb
     mask;
   - the factor assembly F = Q B^QE + T places Q by a plain copy at
     a limb offset;
   - the final division's placement is a limb-offset copy of the
     quotient -- no shift at all.

   An exact bit count is still taken where it is genuinely needed:
   the trimmed slice's mpn_bits feeds the tight series term count.
   The known trade (see the tbound comment in exp_sum_bs.c): the
   splitting integers scale with the frame 64 D rather than the true
   rate r, a content factor (64 D)/r on the leading slices that
   matters for r well below the limb size and fades as r grows;
   deeper slices are limb-aligned in both versions. */
static void
_fixed_exp_reduced_burst(nn_ptr y, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int levels)
{
    slong cap = wn + EXP_BURST_GUARD;
    slong L[FLINT_BITS + 2];
    slong nb = 0, k, nn_, dn;
    slong nexp, dexp, QE;
    nn_ptr num, num2, den, den2, q;
    TMP_INIT;

    /* boundary ladder in limbs: L[0] = max(r/64, 1) (the top slice
       reaches up to the r zero bits regardless), doubling to wn */
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
        /* argument narrower than one boundary step: one slice
           spanning everything */
        L[1] = wn;
        nb = 1;
    }

    {
        slong win_alloc = cap + 1;
        slong q_alloc = wn + 4;

        TMP_START;
        num = TMP_ALLOC((4 * win_alloc + q_alloc) * sizeof(ulong));
        num2 = num + win_alloc;
        den = num2 + win_alloc;
        den2 = den + win_alloc;
        q = den2 + win_alloc;
    }

    num[0] = 1;
    nn_ = 1;
    nexp = 0;
    den[0] = 1;
    dn = 1;
    dexp = 0;
    QE = 0;

    for (k = nb - 1; k >= 0; k--)
    {
        int series = (levels > 0 && k >= levels);
        slong fn, fexp;
        nn_srcptr f;
        TMP_INIT;

        TMP_START;

        if (series)
        {
            /* residual below limb depth L[k]: whole-limb clear */
            nn_ptr xres, fac;

            fac = TMP_ALLOC((wn + 2 + wn) * sizeof(ulong));
            xres = fac + wn + 2;
            flint_mpn_copyi(xres, t, wn);
            flint_mpn_zero(xres + wn - L[k], L[k]);
            fixed_exp_reduced(fac, xres, wn,
                (flint_bitcnt_t) (FLINT_BITS * L[k]), 0);

            f = fac;
            fn = wn + 1;
            fexp = -wn;
        }
        else
        {
            /* slice k = limbs [L[k], L[k+1]) of depth (the top
               slice reaches the top): a pointer and a length */
            slong D = L[k + 1];
            nn_srcptr u = t + (wn - D);
            slong xn = D - (k ? L[k] : 0);
            slong N, tn, qn, QEk;
            nn_ptr T, Q, F;

            while (xn > 0 && u[xn - 1] == 0)
                xn--;
            if (xn == 0)
            {
                TMP_END;
                continue;
            }
            /* strip trailing zero limbs of the slice into the
               frame: u B^-D = (u / B^z) B^-(D - z).  For sparse
               arguments (a low-precision value at high precision:
               one significant limb atop a deep frame) this shrinks
               EVERYTHING downstream -- the kernel's splitting
               integers all scale with D. */
            while (xn > 1 && u[0] == 0)
            {
                u++;
                xn--;
                D--;
            }

            N = _fixed_exp_bs_num_terms((flint_bitcnt_t)
                (FLINT_BITS * D - mpn_bits(u, xn)),
                FLINT_BITS * wn + 64);

            /* kernel buffers per _fixed_exp_sum_bs_powtab_b:
               T spans N (D + 2) + 4 limbs, Q a product of N
               factors below N + 1, F the cap + 3 factor window */
            {
                slong t_alloc = N * (D + 2) + 4;
                slong q_alloc2 = (N * FLINT_BIT_COUNT((ulong) N + 1))
                    / FLINT_BITS + 3;
                slong f_alloc = cap + 3;

                T = TMP_ALLOC((t_alloc + q_alloc2 + f_alloc)
                    * sizeof(ulong));
                Q = T + t_alloc;
                F = Q + q_alloc2;
            }

            _fixed_exp_sum_bs_powtab(T, &tn, Q, &qn, &QEk,
                u, xn, D, N);

            /* DEN *= Q_k, windowed product + swap */
            {
                slong tot = dn + qn;
                slong lo = FLINT_MAX(0, tot - cap);
                flint_mpn_mulmid(den2, den, dn, Q, qn, lo, tot);
                dexp += lo;
                dn = tot - lo;
                while (dn > 1 && den2[dn - 1] == 0)
                    dn--;
                FLINT_SWAP(nn_ptr, den, den2);
            }
            QE += QEk;

            /* factor F = Q B^QEk + T: pure limb placement, exactly
               when it fits within cap + 2 limbs, else only the top
               window with T's tail above the limb wbot (dropping
               the rest loses under one ulp there) */
            {
                slong fl = QEk + qn + 1;    /* F < 2 Q B^QEk */

                FLINT_ASSERT(qn <= cap);
                flint_mpn_zero(F, cap + 3);
                if (fl <= cap + 2)
                {
                    flint_mpn_copyi(F + QEk, Q, qn);
                    if (tn > 0)
                        mpn_add(F, F, cap + 3, T, tn);
                    fexp = 0;
                }
                else
                {
                    slong wbot = fl - (cap + 1);
                    flint_mpn_copyi(F + QEk - wbot, Q, qn);
                    if (wbot < tn)
                    {
                        ulong cy = mpn_add(F, F, cap + 3,
                            T + wbot, tn - wbot);
                        FLINT_ASSERT(cy == 0);
                        (void) cy;
                    }
                    fexp = wbot;
                }
                fn = cap + 3;
                while (fn > 1 && F[fn - 1] == 0)
                    fn--;
                f = F;
            }
        }

        /* NUM *= f, windowed product + swap */
        nexp += fexp;
        {
            slong tot = nn_ + fn;
            slong lo = FLINT_MAX(0, tot - cap);
            flint_mpn_mulmid(num2, num, nn_, f, fn, lo, tot);
            nexp += lo;
            nn_ = tot - lo;
            while (nn_ > 1 && num2[nn_ - 1] == 0)
                nn_--;
            FLINT_SWAP(nn_ptr, num, num2);
        }

        TMP_END;
    }

    /* one balanced Newton division; the placement of the quotient
       is a limb-offset copy, no shift */
    {
        slong E = (nexp + nn_) - (dexp + QE + dn);

        fixed_div_newton(q, num, nn_, den, dn, wn + 2);

        /* y = (q B^-(wn+2)) B^E into wn fraction limbs + units:
           y = q shifted down by 2 - E whole limbs */
        FLINT_ASSERT(2 - E >= 0 && 2 - E <= 3);
        flint_mpn_copyi(y, q + (2 - E), wn + 1);
    }
    TMP_END;
}


void
fixed_exp_reduced(nn_ptr y, nn_srcptr t, slong wn, flint_bitcnt_t r,
    int alg)
{
    FLINT_ASSERT(r >= 16);

    if (alg == 0)
    {
        /* tuned automatic choice; thresholds from tune-exp-reduced.
           For 16 <= r < 32 only the burst paths apply (the series
           functions require t < 2^-32). */
        if (FLINT_BITS * (ulong) wn >= FIXED_EXP_FULLBURST_TERMS * r)
            alg = 4;
        else if (r < 32 || EXP_USE_BURST(wn, r))
            alg = 3;
        else if (EXP_USE_SINH(wn, r))
            alg = 2;
        else
            alg = 1;
    }

    FLINT_ASSERT(alg >= 3 || r >= 32);

    if (alg == 3 || alg == 4)
    {
        _fixed_exp_reduced_burst(y, t, wn, r, (alg == 4) ? 0 : 1);
        return;
    }

    if (alg == 1)
    {
        fixed_exp_rs(y, t, wn);
    }
    else
    {
        nn_ptr s, u2, rt;
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(((wn + 1) + (2 * wn + 1) + (wn + 1))
            * sizeof(ulong));
        u2 = s + (wn + 1);
        rt = u2 + (2 * wn + 1);

        fixed_sinh_rs(s, t, wn);

        /* u2 = (1 + sinh(t)^2) 2^(2 FLINT_BITS wn) */
        flint_mpn_zero(u2, wn);
        flint_mpn_sqrhigh(u2 + wn, s, wn);
        u2[2 * wn] = 1;

        /* cosh(t) = sqrt(1 + sinh(t)^2).  Above the cutoff, run the
           Karp-Markstein square root instead of mpn_sqrtrem: view u2
           as a fraction with 2wn + 2 limbs -- a zero limb on top, the
           widened input normalization of fixed_rsqrt_newton, so that
           the square root's scale factor is a limb shift -- in
           [B^-2, 1.3 B^-2), and compute sqrt with wn + 3 fraction
           limbs.  Truncating the guard limbs leaves rt within
           1 + eps ulps of the floored square root (the two-sided
           error 4 B^(-wn-3) / sqrt(u2') plus half an ulp from the
           short operand stay under one rt ulp) -- absorbed by the
           ~20 ulp reconstruction budget accounted in
           FIXED_EXP_BITWISE_RS_MAX_ERR.

           The default cutoff is deliberately high: up to n = 8192
           this measured as a wash to a few percent slower in context
           on the development machine (mpn_sqrtrem's
           divide-and-conquer being too good and the square root too
           small a fraction of exp; lowering the EXP_USE_SINH
           crossover to exploit a cheaper square root also gained
           nothing).  The Newton square root's asymptotic edge should
           eventually tell; the path is kept behind the cutoff for
           experiments beyond that range.  See dev/notes. */
        if (wn < FIXED_EXP_SQRT_NEWTON_CUTOFF)
        {
            mpn_sqrtrem(rt, NULL, u2, 2 * wn + 1);
        }
        else
        {
            nn_ptr q, a2;
            TMP_INIT;
            TMP_START;
            a2 = TMP_ALLOC((wn + 2 + wn + 5) * sizeof(ulong));
            q = a2 + wn + 2;
            /* the Newton square root reads short input: only the top
               wn + 1 limbs of u2 matter at this precision, plus the
               zero top limb that makes the scale a limb shift */
            flint_mpn_copyi(a2, u2 + wn, wn + 1);
            a2[wn + 1] = 0;
            fixed_sqrt_newton(q, a2, wn + 2, wn + 3);
            flint_mpn_copyi(rt, q + 2, wn + 1);
            TMP_END;
        }

        /* exp(t) = sinh(t) + cosh(t) */
        mpn_add_n(y, rt, s, wn + 1);

        TMP_END;
    }
}
