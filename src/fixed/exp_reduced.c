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
void
_fixed_exp_burst_factor(nn_ptr F, slong * fn, slong * fexp,
    nn_srcptr T, slong tn, nn_srcptr Q, slong qn,
    flint_bitcnt_t Qexp, slong cap)
{
    slong qb = FLINT_BITS * (qn - 1) + FLINT_BIT_COUNT(Q[qn - 1]);
    slong fbits = (slong) Qexp + qb + 1;    /* F < 1.5 Q 2^Qexp */
    slong l;

    /* Q must fit inside the window (in the burst, qbits stays well
       under 0.6 wp against a window of 64 (wn + 4) bits) */
    FLINT_ASSERT(qb < FLINT_BITS * (cap + 1));

    if (fbits <= FLINT_BITS * (cap + 2))
    {
        /* small factor: form it exactly */
        slong Fq = (slong) Qexp / FLINT_BITS;
        int Fb = (int) ((slong) Qexp % FLINT_BITS);

        flint_mpn_zero(F, cap + 3);
        if (Fb)
            F[Fq + qn] = mpn_lshift(F + Fq, Q, qn, Fb);
        else
            flint_mpn_copyi(F + Fq, Q, qn);
        if (tn > 0)
            mpn_add(F, F, cap + 3, T, tn);
        *fexp = 0;
    }
    else
    {
        /* assemble only the top window: wbot is the window's
           bottom bit, limb-aligned */
        slong wbot = ((fbits - FLINT_BITS * (cap + 1))
            / FLINT_BITS) * FLINT_BITS;
        slong sh1 = (slong) Qexp - wbot;
        slong s1q = sh1 / FLINT_BITS;
        int s1b = (int) (sh1 % FLINT_BITS);

        FLINT_ASSERT(sh1 >= 0);
        flint_mpn_zero(F, cap + 3);
        if (s1b)
            F[s1q + qn] = mpn_lshift(F + s1q, Q, qn, s1b);
        else
            flint_mpn_copyi(F + s1q, Q, qn);
        /* add T's limb-aligned tail above wbot; dropping the rest
           loses under one bit at wbot, a relative 2^(-64 cap) */
        {
            slong dq2 = wbot / FLINT_BITS;
            if (dq2 < tn)
            {
                ulong cy = mpn_add(F, F, cap + 3, T + dq2,
                    tn - dq2);
                FLINT_ASSERT(cy == 0);
                (void) cy;
            }
        }
        *fexp = wbot;
    }
    l = cap + 3;
    while (l > 1 && F[l - 1] == 0)
        l--;
    *fn = l;
}

/* Bit-burst evaluation.  Slice boundaries double from depth r:
   b_0 = r, b_{k+1} = min(2 b_k, wp); slice k = bits [b_k, b_{k+1})
   of t, extracted at bit granularity as u = t >> (wp - b_{k+1})
   with the bits above depth b_k masked off (t is read in place; arb
   instead subtracts used bits from a working copy).  Each factor is
   the exact rational

       1 + s_k = (Q_k 2^{Qexp_k} + T_k) / (Q_k 2^{Qexp_k}),

   and instead of dividing per level, the loop accumulates

       NUM = prod (Q_k 2^{Qexp_k} + T_k),      DEN = prod Q_k,
       QE  = sum Qexp_k,

   NUM and DEN msb-truncated to wn + 3 limbs with dropped-bit
   exponents, and finishes with a SINGLE balanced division by
   fixed_div_newton.  The slices are processed in REVERSE (deepest
   first): a deep factor's numerator carries only
   qbits_k + (wp - b_k) + O(1) significant bits -- it is
   Q_k (1 + s_k) with s_k < 2^{-b_k} -- and its denominator only
   qbits_k, so both products grow from small to large and every
   accumulation multiplication is balanced against the content
   gathered so far, the same growing-product scheme as the exp and
   trig reconstructions.  For the level-capped variant (alg 3) the
   remainder below the last splitting slice is evaluated by the
   tuned series choice and folded in as the deepest factor
   (numerator 1 + f_res, denominator 1).

   Frames: every quantity is a mantissa with an explicit power of
   two, value = m 2^e.  A factor F_k = Q_k 2^{Qexp_k} + T_k is
   formed exactly when it fits in the cap window and otherwise only
   its top window is assembled, Q_k shifted into place and T_k's
   limb-aligned tail added, the dropped low bits going into e.
   Errors: each factor and product truncation is one ulp at the
   cap, a relative 2^(-64 (wn + 3) + 64); over at most log2(wp)
   levels on both sides of the quotient, plus the Newton division's
   4 B^{-wn-2}/den and the final one-ulp placement, everything
   lands orders of magnitude inside FIXED_EXP_REDUCED_MAX_ERR.  No
   per-level quotient frames, shifts or additions remain, and the
   product starts from the deepest factor rather than from a
   full-precision 1. */
#define EXP_BURST_GUARD 3

static void
_fixed_exp_reduced_burst(nn_ptr y, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int levels)
{
    slong wp = FLINT_BITS * wn;
    slong cap = wn + EXP_BURST_GUARD;
    slong bnds[FLINT_BITS + 2];
    slong nb = 0, k, nn_, dn, l;
    slong nexp, dexp, QE;
    nn_ptr num, den, sc, q;
    TMP_INIT;

    /* boundary ladder b_0 = r < b_1 < ... < b_nb = wp */
    bnds[nb++] = (slong) r;
    while (bnds[nb - 1] < wp)
    {
        slong nxt = FLINT_MIN(2 * bnds[nb - 1], wp);
        if (levels > 0 && nb > levels)
            nxt = wp;   /* level cap: one final series slice */
        bnds[nb] = nxt;
        nb++;
    }
    nb--;               /* nb slices: [bnds[k], bnds[k+1]) */

    TMP_START;
    num = TMP_ALLOC((2 * (cap + 1) + 2 * cap + 2 + (wn + 4))
        * sizeof(ulong));
    den = num + cap + 1;
    sc = den + cap + 1;
    q = sc + 2 * cap + 2;

    num[0] = 1;
    nn_ = 1;
    nexp = 0;
    den[0] = 1;
    dn = 1;
    dexp = 0;
    QE = 0;

    for (k = nb - 1; k >= 0; k--)
    {
        slong blo = bnds[k], bhi = bnds[k + 1];
        int series = (levels > 0 && k >= levels);
        slong fn, fexp;
        nn_ptr f;
        TMP_INIT;

        TMP_START;

        if (series)
        {
            /* remainder below depth blo: tuned series on a copy
               with the bits above depth blo cleared; factor
               1 + f_res over denominator 1 */
            nn_ptr xres, fac;
            slong zq = blo / FLINT_BITS;
            int zb = (int) (blo % FLINT_BITS);

            fac = TMP_ALLOC((wn + 2 + wn) * sizeof(ulong));
            xres = fac + wn + 2;
            flint_mpn_copyi(xres, t, wn);
            flint_mpn_zero(xres + wn - zq, zq);
            if (zb && zq < wn)
                xres[wn - zq - 1] &= (UWORD(1)
                    << (FLINT_BITS - zb)) - 1;
            fixed_exp_reduced(fac, xres, wn, (flint_bitcnt_t) blo, 0);

            /* value = FacInt 2^(-64 wn), FacInt over wn + 1 limbs
               with the units limb on top; cap >= wn + 1, so no
               truncation is ever needed here */
            f = fac;
            fn = wn + 1;
            fexp = -FLINT_BITS * wn;
        }
        else
        {
            slong sq = (wp - bhi) / FLINT_BITS;
            int sb = (int) ((wp - bhi) % FLINT_BITS);
            slong un = (bhi + FLINT_BITS - 1) / FLINT_BITS + 1;
            slong xn, N, tn, qn, mq;
            int mb;
            flint_bitcnt_t Qexp;
            nn_ptr u, T, Q, F;

            u = TMP_ALLOC(un * sizeof(ulong));
            if (sb)
                mpn_rshift(u, t + sq, wn - sq, sb);
            else
                flint_mpn_copyi(u, t + sq, wn - sq);
            xn = FLINT_MIN(un, wn - sq);
            /* mask the bits above depth blo (shallower slices) */
            mq = (bhi - blo) / FLINT_BITS;
            mb = (int) ((bhi - blo) % FLINT_BITS);
            if (mq < xn)
            {
                if (mb)
                {
                    u[mq] &= (UWORD(1) << mb) - 1;
                    flint_mpn_zero(u + mq + 1, xn - mq - 1);
                }
                else
                    flint_mpn_zero(u + mq, xn - mq);
                xn = FLINT_MIN(xn, mq + 1);
            }
            while (xn > 0 && u[xn - 1] == 0)
                xn--;

            if (xn == 0)
            {
                TMP_END;
                continue;
            }

            {
                slong ubits = FLINT_BITS * (xn - 1)
                    + FLINT_BIT_COUNT(u[xn - 1]);
                N = _fixed_exp_bs_num_terms(
                    (flint_bitcnt_t) (bhi - ubits), wp + 64);
            }

            T = TMP_ALLOC((
                (N * (bhi + 128)) / FLINT_BITS + 4
                + (N * FLINT_BIT_COUNT((ulong) N + 1)) / FLINT_BITS
                + 3 + (cap + 3)) * sizeof(ulong));
            Q = T + (N * (bhi + 128)) / FLINT_BITS + 4;
            F = Q + (N * FLINT_BIT_COUNT((ulong) N + 1))
                / FLINT_BITS + 3;

            _fixed_exp_sum_bs_powtab(T, &tn, Q, &qn, &Qexp,
                u, xn, (flint_bitcnt_t) bhi, N);

            /* DEN *= Q_k, msb-truncated at cap */
            if (dn >= qn)
                flint_mpn_mul(sc, den, dn, Q, qn);
            else
                flint_mpn_mul(sc, Q, qn, den, dn);
            l = dn + qn;
            while (l > 1 && sc[l - 1] == 0)
                l--;
            if (l > cap)
            {
                flint_mpn_copyi(den, sc + (l - cap), cap);
                dexp += FLINT_BITS * (l - cap);
                dn = cap;
            }
            else
            {
                flint_mpn_copyi(den, sc, l);
                dn = l;
            }
            QE += (slong) Qexp;

            _fixed_exp_burst_factor(F, &fn, &fexp, T, tn, Q, qn,
                Qexp, cap);
            f = F;
        }

        /* NUM *= f, msb-truncated at cap */
        if (nn_ >= fn)
            flint_mpn_mul(sc, num, nn_, f, fn);
        else
            flint_mpn_mul(sc, f, fn, num, nn_);
        l = nn_ + fn;
        while (l > 1 && sc[l - 1] == 0)
            l--;
        nexp += fexp;
        if (l > cap)
        {
            flint_mpn_copyi(num, sc + (l - cap), cap);
            nexp += FLINT_BITS * (l - cap);
            nn_ = cap;
        }
        else
        {
            flint_mpn_copyi(num, sc, l);
            nn_ = l;
        }

        TMP_END;
    }

    /* y = (num 2^nexp) / (den 2^{dexp + QE}), one balanced Newton
       division.  Normalize both mantissas to fractions in [1/2, 1)
       (so the denominator's top limb is nonzero as the kernel
       requires); the log difference is then exactly
       (nexp + numb) - (dexp + QE + denb) with numb, denb the
       mantissa bit lengths, and y in [1, e) makes that difference
       0 or 1. */
    {
        slong numb = FLINT_BITS * (nn_ - 1)
            + FLINT_BIT_COUNT(num[nn_ - 1]);
        slong denb = FLINT_BITS * (dn - 1)
            + FLINT_BIT_COUNT(den[dn - 1]);
        slong sh = (nexp + numb) - (dexp + QE + denb);
        slong tot, yq;
        int lz, yb;
        nn_ptr a2, b2;

        a2 = sc;
        b2 = sc + cap + 1;

        lz = (int) ((FLINT_BITS - (denb % FLINT_BITS))
            % FLINT_BITS);
        if (lz)
            mpn_lshift(a2, den, dn, lz);
        else
            flint_mpn_copyi(a2, den, dn);

        lz = (int) ((FLINT_BITS - (numb % FLINT_BITS))
            % FLINT_BITS);
        if (lz)
            mpn_lshift(b2, num, nn_, lz);
        else
            flint_mpn_copyi(b2, num, nn_);

        fixed_div_newton(q, b2, nn_, a2, dn, wn + 2);

        /* q = b/a in [1/2, 2) over wn + 2 fraction limbs; value
           y = (b/a) 2^sh with sh in {0, 1}.  Place into the y
           frame (wn fraction limbs + units): a plain truncating
           shift, one-sided within 1 ulp on top of the division's
           error. */
        FLINT_ASSERT(sh == 0 || sh == 1);
        tot = 2 * FLINT_BITS - sh;
        yq = tot / FLINT_BITS;
        yb = (int) (tot % FLINT_BITS);
        if (yb)
        {
            mpn_rshift(q, q + yq, wn + 4 - yq, yb);
            flint_mpn_copyi(y, q, wn + 1);
        }
        else
            flint_mpn_copyi(y, q + yq, wn + 1);
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
        nn_ptr s, u2, rt, rem;
        TMP_INIT;

        TMP_START;
        s = TMP_ALLOC(((wn + 1) + (2 * wn + 1) + (wn + 1) + (wn + 2))
            * sizeof(ulong));
        u2 = s + (wn + 1);
        rt = u2 + (2 * wn + 1);
        rem = rt + (wn + 1);

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
