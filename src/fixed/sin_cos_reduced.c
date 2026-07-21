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

   The burst machinery mirrors _fixed_exp_reduced_burst at LIMB
   granularity throughout: the boundary ladder is a list of limb
   counts (tripling), a slice is a pointer + length into t with a
   whole-limb mask, dead low limbs of a slice fold into its frame
   (u B^-D = (u/B^z) B^-(D-z)), and each slice's factor is the
   exact Gaussian rational

       exp(i x_k) ~ (Q_k B^{QE_k} - A_k + i B_k) / (Q_k B^{QE_k})

   truncated at N_k series terms, all exponents limb counts.
   Instead of combining per level, the loop accumulates the complex
   numerator NUM = prod (Q_k B^{QE_k} - A_k + i B_k) and the real
   denominator DEN = prod Q_k, both as WINDOWED middle products
   msb-truncated to wn + 3 limbs with dropped-limb exponents
   (mulmid + ping-pong for DEN; the two-regime complex update in
   trig_num_mul for NUM), deepest slice first so every complex
   accumulation multiplication is balanced against the content
   gathered so far.  The finish is TWO balanced Newton divisions
   (sine and cosine) against the single accumulated denominator,
   fed unnormalized (fixed_div_newton needs only a nonzero top
   limb) with quotient placement a limb-offset copy -- where
   arb_sin_cos_arf_bb spends one full-precision square root PER
   SLICE (cos_k from sin_k) plus a full-precision complex product
   tree.  g comes from the cosine quotient by negation: the ~2r-bit
   cancellation is harmless because g is only needed to absolute
   2^(-64 wn), like the direct series path.

   Errors: each factor window and product drop is one-sided ulps at
   the cap per component, with the mulmid boundary slack held three
   limbs below the kept frames; over at most log2(wn) + 1 levels,
   plus the Newton divisions' 4 B^{-wn-2}/den and the final
   one-ulp placements, everything lands far inside
   FIXED_SIN_COS_REDUCED_MAX_ERR = 96 ulps, the same budget as
   fixed_exp_reduced. */

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
   development VM. */
#ifndef FIXED_TRIG_REDUCED_SINSQRT_TERMS
#define FIXED_TRIG_REDUCED_SINSQRT_TERMS 50
#endif
#ifndef FIXED_TRIG_BURST_TERMS_SMALL_R
#define FIXED_TRIG_BURST_TERMS_SMALL_R 320
#endif
#ifndef FIXED_TRIG_BURST_TERMS
#define FIXED_TRIG_BURST_TERMS 2560
#endif
#ifndef FIXED_TRIG_FULLBURST_TERMS
#define FIXED_TRIG_FULLBURST_TERMS 393216
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
    (FLINT_BITS * (wn) >= FIXED_TRIG_REDUCED_SINSQRT_TERMS * (slong) (r))
#define TRIG_USE_BURST(wn, r) \
    (FLINT_BITS * (wn) >= (((r) < 64) ? FIXED_TRIG_BURST_TERMS_SMALL_R \
        : FIXED_TRIG_BURST_TERMS) * (slong) (r))

#ifndef FIXED_TRIG_REDUCED_SQRT_NEWTON_CUTOFF
#define FIXED_TRIG_REDUCED_SQRT_NEWTON_CUTOFF 2000
#endif

/* g = 1 - cos = 1 - sqrt(1 - sin^2) from the wn-limb sine
   fraction ss: a squaring and a square root, the reconstruction
   shared by algorithm 2 and the burst's sine-only slices.  The
   square root sits just below 1, so its derivative is ~1/2 and
   nothing is amplified: sqrhigh's few-ulp deficit on sin^2, the
   square root's floor and the Newton route's ~2 compensated ulps
   leave g within a few 2^-64wn ulps either way. */
static void
_fixed_g_from_sin(nn_ptr yg, nn_srcptr ss, slong wn)
{
    nn_ptr u2, rt;
    TMP_INIT;

    TMP_START;
    u2 = TMP_ALLOC(((2 * wn + 2) + (wn + 4)) * sizeof(ulong));
    rt = u2 + (2 * wn + 2);

    flint_mpn_sqrhigh(rt, ss, wn);
    if (flint_mpn_zero_p(rt, wn))
    {
        flint_mpn_zero(yg, wn);         /* sin^2 below 2^-64wn */
    }
    else if (wn < FIXED_TRIG_REDUCED_SQRT_NEWTON_CUTOFF)
    {
        flint_mpn_zero(u2, wn);
        mpn_neg(u2 + wn, rt, wn);       /* (1 - sin^2) 2^(128wn) */
        mpn_sqrtrem(rt, NULL, u2, 2 * wn);
        mpn_neg(yg, rt, wn);            /* g 2^(64wn) */
    }
    else
    {
        /* the Newton square root reads short input directly:
           (1 - sin^2) as a wn-limb fraction, no zero padding */
        mpn_neg(u2, rt, wn);
        fixed_sqrt_newton(rt, u2, wn, wn + 2);
        if (rt[wn + 2])
            flint_mpn_zero(yg, wn);     /* cos rounded to 1 */
        else
            mpn_neg(yg, rt + 2, wn);
    }
    TMP_END;
}

static slong
nnn_add_local(nn_ptr z, slong zn, nn_srcptr a, slong an)
{
    ulong cy;
    if (an == 0)
        return zn;
    if (zn >= an)
    {
        cy = mpn_add(z, z, zn, a, an);
        z[zn] = cy;
        return zn + (cy != 0);
    }
    flint_mpn_zero(z + zn, an - zn);
    cy = mpn_add(z, a, an, z, zn);
    z[an] = cy;
    return an + (cy != 0);
}

/* Window bottom (in LIMBS) for a slice factor whose real part is
   R = Q B^QE - A ~ Q B^QE: 0 when the factor fits within cap + 2
   limbs exactly, else so the kept window spans about cap + 1
   limbs. */
static slong
_trig_factor_wbot(slong qn, slong QE, slong cap)
{
    slong fl = QE + qn;

    FLINT_ASSERT(qn <= cap + 1);
    if (fl <= cap + 2)
        return 0;
    return fl - (cap + 1);
}

/* F = Q B^(QE - wbot) - U B^(uexp - wbot): the common assembly of
   both components of the slice factor (cosine with U = A, the sine
   prefactor with U = B), all frames in LIMBS.  Q is placed by a
   plain copy at its offset and U is subtracted AT ITS OFFSET
   directly (dominance-guaranteed) -- no scratch, no shifts.  When
   U's frame sits below wbot its low limbs are dropped, understating
   the subtraction by under one ulp at wbot, a relative
   2^(-64 cap), like the exp factor's dropped tail.  F must have
   room for cap + 3 limbs. */
static void
_trig_qminus(nn_ptr F, slong * fn, nn_srcptr U, slong un,
    slong uexp, nn_srcptr Q, slong qn, slong QE, slong wbot,
    slong cap)
{
    slong sh = uexp - wbot;
    slong l;

    FLINT_ASSERT(QE >= wbot);
    flint_mpn_zero(F, cap + 3);
    flint_mpn_copyi(F + (QE - wbot), Q, qn);

    if (sh >= 0)
    {
        FLINT_ASSERT(un + sh <= cap + 3);
        if (un > 0)
        {
            ulong bw = mpn_sub(F + sh, F + sh, cap + 3 - sh,
                U, un);
            FLINT_ASSERT(bw == 0);
            (void) bw;
        }
    }
    else if (un + sh > 0)
    {
        /* drop U's low (-sh) limbs */
        ulong bw = mpn_sub(F, F, cap + 3, U - sh, un + sh);
        FLINT_ASSERT(bw == 0);
        (void) bw;
    }

    l = cap + 3;
    while (l > 1 && F[l - 1] == 0)
        l--;
    *fn = l;
}

static void
trig_num_mul(nn_ptr nc, slong * ncn, nn_ptr ns, slong * nsn,
    slong * nexp, nn_srcptr fc, slong fcn, nn_srcptr fs, slong fsn,
    slong fexp, nn_ptr sc, slong cap)
{
    /* Complex update in mpn magnitudes over WINDOWED middle
       products with a common limb bottom lo -- only the surviving
       limbs of each product are computed, the mulmid boundary
       slack (a deficit under ~B^(lo+3)) and the drop being
       one-sided ulps at lo, the same class as the factor windows.

       Two regimes, keyed to the size of the sine-sine product:

       - Karatsuba (arb's pmerge) when tot2 = nsn + fsn >= lo + 8:
         m1 = nc fc, m2 = ns fs, m3 = (nc + ns)(fc + fs);
         re = m1 - m2, im = m3 - m1 - m2.  The windowed
         subtractions need the true differences to dominate the
         slack: im >= ns fc >= B^(lo + 6) under the threshold
         (fcn >= fsn since cos > sin), three limbs above it, and
         re is always of full magnitude.

       - direct otherwise (short sine operands, so the extra
         products are cheap): re = m1 - m2 with m2 computed in
         full (it is tiny), and im = nc fs + ns fc + ns fs
         accumulated as a SUM of nonnegative windows -- additions
         cannot borrow, which is what makes the near-cancelling
         regime safe.

       All frames are limb counts; both components end
       msb-truncated at cap in a common frame keyed to the real
       part, which always dominates (|arg| < 2^-15). */
    nn_ptr m1 = sc, m2 = sc + (cap + 8),
           m3 = sc + 2 * (cap + 8),
           s1 = sc + 3 * (cap + 8),
           s2 = s1 + (cap + 8);
    slong tot1 = *ncn + fcn, tot2 = *nsn + fsn, tot3, lo;
    slong l1, l2, l3, ls1, ls2, lmx, drop;
    ulong bw;

    /* three limbs of margin below the natural bottom: the mulmid
       boundary slack then sits entirely below the frame the final
       truncation keeps, so it stays one-sided ulps there even for
       the imaginary component, whose own top can run a few limbs
       short of the real part's */
    lo = FLINT_MAX(0, tot1 - cap - 3);

    flint_mpn_mulmid(m1, nc, *ncn, fc, fcn, lo, tot1);
    l1 = tot1 - lo;
    while (l1 > 1 && m1[l1 - 1] == 0)
        l1--;

    if (*nsn > 0 && fsn > 0 && tot2 >= lo + 8)
    {
        /* Karatsuba: s1 = nc + ns, s2 = fc + fs,
           m3 = s1 s2, m2 = ns fs */
        flint_mpn_copyi(s1, nc, *ncn);
        ls1 = nnn_add_local(s1, *ncn, ns, *nsn);
        flint_mpn_copyi(s2, fc, fcn);
        ls2 = nnn_add_local(s2, fcn, fs, fsn);
        tot3 = ls1 + ls2;
        FLINT_ASSERT(tot3 <= tot1 + 2);

        flint_mpn_mulmid(m2, ns, *nsn, fs, fsn, lo, tot2);
        l2 = tot2 - lo;
        while (l2 > 1 && m2[l2 - 1] == 0)
            l2--;

        flint_mpn_mulmid(m3, s1, ls1, s2, ls2, lo, tot3);
        l3 = tot3 - lo;
        while (l3 > 1 && m3[l3 - 1] == 0)
            l3--;

        /* im: m3 -= m1 + m2; re: m1 -= m2 */
        bw = mpn_sub(m3, m3, l3, m1, l1);
        FLINT_ASSERT(bw == 0);
        bw = mpn_sub(m3, m3, l3, m2, l2);
        FLINT_ASSERT(bw == 0);
        bw = mpn_sub(m1, m1, l1, m2, l2);
        FLINT_ASSERT(bw == 0);
        (void) bw;
        while (l3 > 1 && m3[l3 - 1] == 0)
            l3--;
        while (l1 > 1 && m1[l1 - 1] == 0)
            l1--;
    }
    else
    {
        /* direct: im = nc fs + ns fc as nonnegative windows at
           the common bottom (additions cannot borrow), and
           re -= ns fs with its tiny windowed product */
        slong tota = *ncn + fsn, totb = *nsn + fcn;

        l3 = 0;
        if (fsn > 0 && tota > lo)
        {
            flint_mpn_mulmid(m3, nc, *ncn, fs, fsn, lo, tota);
            l3 = tota - lo;
            while (l3 > 1 && m3[l3 - 1] == 0)
                l3--;
        }
        if (*nsn > 0 && totb > lo)
        {
            flint_mpn_mulmid(s1, ns, *nsn, fc, fcn, lo, totb);
            ls1 = totb - lo;
            while (ls1 > 1 && s1[ls1 - 1] == 0)
                ls1--;
            if (l3 == 0)
            {
                flint_mpn_copyi(m3, s1, ls1);
                l3 = ls1;
            }
            else
                l3 = nnn_add_local(m3, l3, s1, ls1);
        }
        if (*nsn > 0 && fsn > 0 && tot2 > lo)
        {
            /* the sine-sine term: windowed too (l2 <= 7 under the
               regime threshold); its mulmid deficit understates
               both the im addition and the re subtraction --
               one-sided, the right way */
            flint_mpn_mulmid(s2, ns, *nsn, fs, fsn, lo, tot2);
            l2 = tot2 - lo;
            while (l2 > 1 && s2[l2 - 1] == 0)
                l2--;
            if (!(l2 == 1 && s2[0] == 0))
            {
                bw = mpn_sub(m1, m1, l1, s2, l2);
                FLINT_ASSERT(bw == 0);
                (void) bw;
                while (l1 > 1 && m1[l1 - 1] == 0)
                    l1--;
            }
        }
    }

    /* common truncation keyed to the real part, frames in LIMBS */
    *nexp += fexp + lo;
    lmx = FLINT_MAX(l1, l3);
    FLINT_ASSERT(lmx == l1 || l3 <= l1 + 1);
    drop = (lmx > cap) ? (lmx - cap) : 0;
    if (drop > 0)
    {
        flint_mpn_copyi(nc, m1 + drop, FLINT_MAX(1, l1 - drop));
        *ncn = FLINT_MAX(1, l1 - drop);
        if (l3 > drop)
        {
            flint_mpn_copyi(ns, m3 + drop, l3 - drop);
            *nsn = l3 - drop;
        }
        else
        {
            ns[0] = 0;
            *nsn = 0;
        }
        *nexp += drop;
    }
    else
    {
        flint_mpn_copyi(nc, m1, l1);
        *ncn = l1;
        if (l3 > 0)
            flint_mpn_copyi(ns, m3, l3);
        *nsn = l3;
    }
}

/* one balanced Newton division NUM B^nexp / (DEN B^dexp), placed
   into a wn-limb (unit == 0) or (wn+1)-limb frame; mirrors the
   tail of _fixed_exp_reduced_burst: fixed_div_newton needs only a
   nonzero top LIMB on the denominator, so neither operand is
   normalized -- the quotient's frame is recovered from the limb
   counts alone, E = (nexp + nn_) - (dexp + dn) in [-2, 1]-ish for
   the cosine and anywhere below for a tiny sine, and the output
   placement is a limb copy.  y must have wn + 1 limbs of room;
   q needs wn + 4. */
static void
trig_burst_div(nn_ptr y, nn_srcptr num, slong nn_, slong nexp,
    nn_srcptr den, slong dn, slong dexp, slong wn, nn_ptr q)
{
    slong E, yq, avail;

    if (nn_ == 0 || (nn_ == 1 && num[0] == 0))
    {
        flint_mpn_zero(y, wn + 1);
        return;
    }

    FLINT_ASSERT(den[dn - 1] != 0);

    fixed_div_newton(q, num, nn_, den, dn, wn + 2);

    /* q = num/den in (0, B) over wn + 2 fraction limbs + a units
       limb at q[wn + 2]; value y = (num/den) B^E.  sin, cos <= 1,
       so E <= 1 always; a tiny sine can push E far negative.
       Place into wn fraction limbs + a units limb by dropping the
       bottom (2 - E) limbs of the quotient frame. */
    E = (nexp + nn_) - (dexp + dn);
    FLINT_ASSERT(E <= 1);
    yq = 2 - E;
    if (yq > wn + 3)
    {
        /* value below one ulp of the output frame */
        flint_mpn_zero(y, wn + 1);
        return;
    }
    avail = wn + 4 - yq;        /* limbs of q above the drop */
    flint_mpn_copyi(y, q + yq, FLINT_MIN(wn + 1, avail));
    if (avail < wn + 1)
        flint_mpn_zero(y + avail, wn + 1 - avail);
}

static void
_fixed_sin_cos_reduced_burst(nn_ptr ysin, nn_ptr yg, nn_srcptr t,
    slong wn, flint_bitcnt_t r, int levels)
{
    slong cap = wn + TRIG_BURST_GUARD;
    slong L[FLINT_BITS + 2];
    slong nb = 0, k, ncn, nsn, dn;
    slong nexp, dexp, QE;
    nn_ptr nc, ns, den, den2, sc, q, ycos;
    TMP_INIT;

    /* boundary ladder in LIMBS, TRIPLING (like
       arb_sin_cos_arf_bb's bits *= 3): against doubling this cuts
       the level count by lg 3 and the total splitting-tree content
       from 2x to 1.5x the dominant first tree, at slightly larger
       per-slice trees.  L[0] = max(r/64, 1): the top slice reaches
       up to the r zero bits regardless. */
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
        /* argument narrower than one boundary step: one slice
           spanning everything */
        L[1] = wn;
        nb = 1;
    }

    TMP_START;
    nc = TMP_ALLOC((4 * (cap + 1)
        + 5 * (cap + 8) + (wn + 4)
        + (wn + 1)) * sizeof(ulong));
    ns = nc + (cap + 1);
    den = ns + (cap + 1);
    den2 = den + (cap + 1);
    sc = den2 + (cap + 1);
    q = sc + 5 * (cap + 8);
    ycos = q + (wn + 4);

    nc[0] = 1;
    ncn = 1;
    ns[0] = 0;
    nsn = 0;
    nexp = 0;
    den[0] = 1;
    dn = 1;
    dexp = 0;
    QE = 0;

    for (k = nb - 1; k >= 0; k--)
    {
        int series = (levels > 0 && k >= levels);
        slong fcn, fsn, fexp;
        nn_srcptr fc, fs;
        TMP_INIT;

        TMP_START;

        if (series)
        {
            /* remainder below depth 64 L[k] by the tuned choice;
               factor (1 - g_res) + i s_res over denominator 1.
               The residual mask is a whole-limb clear. */
            nn_ptr fac, facs, xres;

            fac = TMP_ALLOC((3 * (wn + 2)) * sizeof(ulong));
            facs = fac + (wn + 2);
            xres = facs + (wn + 2);
            flint_mpn_copyi(xres, t, wn);
            flint_mpn_zero(xres + wn - L[k], L[k]);
            fixed_sin_cos_reduced(facs, fac, xres, wn,
                (flint_bitcnt_t) (FLINT_BITS * L[k]), 0);

            /* FC = B^wn - g_res, FS = s_res */
            mpn_neg(fac, fac, wn);
            fac[wn] = flint_mpn_zero_p(fac, wn);
            fcn = wn + 1;
            while (fcn > 1 && fac[fcn - 1] == 0)
                fcn--;
            fc = fac;
            fsn = wn;
            while (fsn > 0 && facs[fsn - 1] == 0)
                fsn--;
            fs = facs;
            fexp = -wn;
        }
        else
        {
            slong D = L[k + 1];
            nn_srcptr u = t + (wn - D);
            slong xn = D - (k ? L[k] : 0);
            slong N, an, bn2, qn, ae, be, QEk;
            nn_ptr A, B, Q, FC, FS;

            while (xn > 0 && u[xn - 1] == 0)
                xn--;
            if (xn == 0)
            {
                TMP_END;
                continue;
            }
            /* strip trailing zero limbs of the slice into the
               frame: u B^-D = (u / B^z) B^-(D - z), shrinking the
               whole tree for sparse arguments */
            while (xn > 1 && u[0] == 0)
            {
                u++;
                xn--;
                D--;
            }

            {
                slong ubits = FLINT_BITS * (xn - 1)
                    + FLINT_BIT_COUNT(u[xn - 1]);
                /* the joint tree runs over the y = x^2 series:
                   half the exp count (the (2j+1)! denominators
                   only shrink terms further), re-padded to a high
                   2-valuation for the fixed midpoint splitting */
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
                A = TMP_ALLOC((2 * (cap + 2 + qb2 + 8) + qb2
                    + 2 * (cap + 5)) * sizeof(ulong));
                B = A + (cap + 2 + qb2 + 8);
                Q = B + (cap + 2 + qb2 + 8);
                FC = Q + qb2;
                FS = FC + (cap + 5);
            }

            {
                int slice_sqrt =
                    (N >= FIXED_TRIG_SLICE_SQRT_TERMS);
                _fixed_sin_cos_sum_bs_powtab(
                    slice_sqrt ? NULL : A, &an, &ae, B, &bn2, &be,
                    Q, &qn, &QEk, u, xn, D, N, cap + 2);
                if (slice_sqrt)
                    an = -1;    /* flag for the factor assembly */
            }

            /* DEN *= Q_k: windowed product + swap (ping-pong) */
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

            /* factor (over the single denominator Q B^QEk):
               FC = Q B^QEk - A            (~ cos_k Q B^QEk)
               FS = x (Q B^QEk - B) B^-D   (= sin_k Q B^QEk)
               assembled in a common window with bottom wbot; the
               sine's B^-D is a truncating limb drop of the small
               product against the slice integer x */
            {
                slong wbot = _trig_factor_wbot(qn, QEk, cap);
                slong w2n, pl;
                nn_ptr prod;
                TMP_INIT;

                if (an >= 0)
                    _trig_qminus(FC, &fcn, A, an, ae, Q, qn, QEk,
                        wbot, cap);
                _trig_qminus(FS, &w2n, B, bn2, be, Q, qn, QEk,
                    wbot, cap);

                TMP_START;
                prod = TMP_ALLOC((w2n + xn + 2) * sizeof(ulong));
                if (w2n >= xn)
                    flint_mpn_mul(prod, FS, w2n, u, xn);
                else
                    flint_mpn_mul(prod, u, xn, FS, w2n);
                pl = w2n + xn;
                while (pl > 1 && prod[pl - 1] == 0)
                    pl--;
                if (D >= pl)
                {
                    FS[0] = 0;
                    fsn = 0;
                }
                else
                {
                    flint_mpn_copyi(FS, prod + D, pl - D);
                    fsn = pl - D;
                    while (fsn > 0 && FS[fsn - 1] == 0)
                        fsn--;
                }
                if (an < 0)
                {
                    /* sine-only slice, cosine window on the SAME
                       denominator: with
                       FS = sin_k Q B^(QEk - wbot) as an integer
                       window and Dq = Q B^sh1 (sh1 = QEk - wbot),

                           FC = cos_k Q B^(QEk - wbot)
                              = sqrt(Dq^2 - FS^2)

                       -- no division at all: Dq^2 is a short exact
                       squaring placed at a limb offset, FS^2 one
                       full squaring of the window, and a single
                       square root finishes.  FS's own window error
                       delta contributes 2 FS delta / (2 FC)
                       ~ delta ulps through the square root, so all
                       errors stay a few one-sided ulps at the
                       window bottom, the same class as the joint
                       tree's truncations. */
                    slong sh1 = QEk - wbot;

                    while (fsn > 0 && FS[fsn - 1] == 0)
                        fsn--;

                    if (fsn == 0)
                    {
                        /* zero sine window: FC = Q B^sh1 */
                        flint_mpn_zero(FC, sh1);
                        flint_mpn_copyi(FC + sh1, Q, qn);
                        fcn = sh1 + qn;
                        while (fcn > 1 && FC[fcn - 1] == 0)
                            fcn--;
                    }
                    else
                    {
                        nn_ptr q2, v, f2, w;
                        slong q2n, vlen, vn, f2n;

                        q2 = TMP_ALLOC(((2 * qn + 2)
                            + (2 * cap + 2 * qn + 12)
                            + (2 * cap + 8) + (2 * cap + 12))
                            * sizeof(ulong));
                        v = q2 + (2 * qn + 2);
                        f2 = v + (2 * cap + 2 * qn + 12);
                        w = f2 + (2 * cap + 8);

                        /* V = Dq^2 = Q^2 B^(2 sh1) */
                        flint_mpn_sqr(q2, Q, qn);
                        q2n = 2 * qn;
                        while (q2n > 1 && q2[q2n - 1] == 0)
                            q2n--;
                        flint_mpn_zero(v, 2 * sh1);
                        flint_mpn_copyi(v + 2 * sh1, q2, q2n);
                        vlen = 2 * sh1 + q2n;

                        /* V -= FS^2 (FS <= Dq, so it fits) */
                        flint_mpn_sqr(f2, FS, fsn);
                        f2n = 2 * fsn;
                        while (f2n > 1 && f2[f2n - 1] == 0)
                            f2n--;
                        FLINT_ASSERT(f2n <= vlen);
                        {
                            ulong bw = mpn_sub(v, v, vlen, f2, f2n);
                            FLINT_ASSERT(bw == 0);
                            (void) bw;
                        }
                        vn = vlen;
                        while (vn > 1 && v[vn - 1] == 0)
                            vn--;

                        if (vn < 2 * FIXED_TRIG_REDUCED_SQRT_NEWTON_CUTOFF)
                        {
                            /* exact integer square root */
                            fcn = (vn + 1) / 2;
                            mpn_sqrtrem(FC, NULL, v, vn);
                            while (fcn > 1 && FC[fcn - 1] == 0)
                                fcn--;
                        }
                        else
                        {
                            /* Newton at pure limb granularity:
                               take V's top limbs at an EVEN limb
                               position -- vh = ceil(vn / 2) * 2 --
                               giving vhat in [B^-2, 1) over nin
                               limbs (a zero top limb when vn is
                               odd), which is exactly
                               fixed_sqrt_newton's accepted range;
                               then sqrt(V) = sqrt(vhat) B^(vh/2)
                               and the result placement is a limb
                               copy at offset vh/2 - (cap + 1). */
                            slong vh = vn + (vn & 1);
                            slong nin = cap + 3;
                            slong off, avail, e;
                            nn_ptr vf = w, rt = w + nin;

                            off = vh - nin;
                            FLINT_ASSERT(off > 0);
                            avail = vn - off;
                            flint_mpn_copyi(vf, v + off, avail);
                            if (avail < nin)
                                flint_mpn_zero(vf + avail,
                                    nin - avail);

                            fixed_sqrt_newton(rt, vf, nin, cap + 1);
                            /* rt: (cap + 1)-limb fraction of
                               sqrt(vhat) with a units limb at
                               rt[cap + 1];
                               FC = rt-as-integer
                                    * B^(vh/2 - (cap + 1)) */
                            e = vh / 2 - (cap + 1);
                            if (e >= 0)
                            {
                                flint_mpn_zero(FC, e);
                                flint_mpn_copyi(FC + e, rt, cap + 2);
                                fcn = e + cap + 2;
                            }
                            else
                            {
                                flint_mpn_copyi(FC, rt + (-e),
                                    cap + 2 - (-e));
                                fcn = cap + 2 - (-e);
                            }
                            while (fcn > 1 && FC[fcn - 1] == 0)
                                fcn--;
                        }
                    }
                }

                TMP_END;
                fexp = wbot;
            }
            fc = FC;
            fs = FS;
        }

        trig_num_mul(nc, &ncn, ns, &nsn, &nexp, fc, fcn, fs, fsn,
            fexp, sc, cap);

        TMP_END;
    }

    /* two balanced Newton divisions against the common denominator */
    trig_burst_div(ysin, ns, nsn, nexp, den, dn, dexp + QE, wn, q);
    trig_burst_div(ycos, nc, ncn, nexp, den, dn, dexp + QE, wn, q);

    /* g = 1 - cos from the cosine quotient */
    if (ycos[wn])
        flint_mpn_zero(yg, wn);
    else
        mpn_neg(yg, ycos, wn);

    /* ysin's units limb was scratch in the division frame */
    TMP_END;
}

void
fixed_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg)
{
    FLINT_ASSERT(r >= 16);

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

    FLINT_ASSERT(alg >= 3 || r >= 32);

    if (alg == 3 || alg == 4)
    {
        nn_ptr ys;
        TMP_INIT;
        TMP_START;
        /* the divisions want a units limb of room */
        ys = TMP_ALLOC((wn + 1) * sizeof(ulong));
        _fixed_sin_cos_reduced_burst(ys, yg, t, wn, r,
            (alg == 4) ? 0 : 1);
        flint_mpn_copyi(ysin, ys, wn);
        TMP_END;
        return;
    }

    if (alg == 2)
    {
        /* sine series + square root, factored out of
           _fixed_tan_halfangle_mid */
        nn_ptr ss;
        TMP_INIT;

        TMP_START;
        ss = TMP_ALLOC((wn + 2) * sizeof(ulong));
        fixed_sin_rs(ss, t, wn);
        flint_mpn_copyi(ysin, ss, wn);
        _fixed_g_from_sin(yg, ss, wn);
        TMP_END;
    }
    else
    {
        /* both series; factored out of _fixed_tan_halfangle_mid */
        nn_ptr ss, cc;
        TMP_INIT;

        TMP_START;
        ss = TMP_ALLOC(2 * (wn + 2) * sizeof(ulong));
        cc = ss + (wn + 2);

        fixed_sin_cos_rs(ss, cc, t, wn);
        flint_mpn_copyi(ysin, ss, wn);

        if (cc[wn])
            flint_mpn_zero(yg, wn);     /* cos t = 1: g = 0 */
        else
            mpn_neg(yg, cc, wn);        /* g = 1 - cos t */
        TMP_END;
    }
}
