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
   (mulmid + ping-pong for DEN; flint_mpn_mulhigh_n_complex in
   trig_num_mul for NUM, which shares transforms across its three
   products above the complex FFT cutoff), deepest slice first so
   every complex accumulation multiplication is balanced against
   the content gathered so far.  The finish is TWO balanced Newton divisions
   (sine and cosine) against the single accumulated denominator,
   fed unnormalized (fixed_div_newton needs only a nonzero top
   limb) with quotient placement a limb-offset copy -- where
   arb_sin_cos_arf_bb spends one full-precision square root PER
   SLICE (cos_k from sin_k) plus a full-precision complex product
   tree.  g comes from the cosine quotient by negation: the ~2r-bit
   cancellation is harmless because g is only needed to absolute
   2^(-64 wn), like the direct series path.

   Errors: each factor window and product drop is a few ulps at the
   cap per component (one-sided for the mulmid windows, whose
   boundary slack is held three limbs below the kept frames;
   two-sided, below 3 ulps of the lowest returned limb, for the high
   complex products of the NUM updates); over at most log2(wn) + 1
   levels,
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
    slong fexp, slong cap)
{
    /* Complex update NUM <- NUM F through flint_mpn_mulhigh_n_complex.
       Both pairs are brought to a common length n = max(ncn, fcn)
       by padding the shorter pair BELOW with zeros (its exponent
       absorbs the pad: the parts of a pair share one exponent, so
       the sine part is padded like its cosine part and zero-extended
       on top), the high product gives limbs [n, 2n] of each part --
       the top n + 1 limbs of a product whose top limb is nonzero
       since both cosine parts are -- and the frame is closed by the
       msb truncation to cap keyed to the real part, which dominates
       (|arg| < 2^-15).  The parts are magnitudes: the real part
       nc fc - ns fs is positive of full size, the imaginary part
       nc fs + ns fc nonnegative; a negative sign can only be the
       two-sided rounding (below 3 ulps of the lowest limb) of a
       vanishing imaginary part and is read as zero.  Above
       flint_mpn_mul_complex_fft_cutoff limbs the three products
       share their transforms. */
    slong n, pn, pf, k, drop, lr, li;
    nn_ptr ar, ai, br, bi, rr, ri;
    int sr, si;
    TMP_INIT;

    /* NUM = 1 (the first factor): NUM F = F exactly -- a high product
       would drop F's lowest limb, one of the frame's three guard
       limbs, before anything is gained */
    if (*ncn == 1 && nc[0] == 1 && (*nsn == 0 || (*nsn == 1 && ns[0] == 0)))
    {
        /* copy F truncated to cap limbs keyed to its real part (the
           factor can be longer than the accumulator's frame) */
        drop = (fcn > cap) ? (fcn - cap) : 0;
        flint_mpn_copyi(nc, fc + drop, fcn - drop);
        *ncn = fcn - drop;
        if (fsn > drop)
        {
            flint_mpn_copyi(ns, fs + drop, fsn - drop);
            *nsn = fsn - drop;
        }
        else
        {
            ns[0] = 0;
            *nsn = 1;
        }
        *nexp += fexp + drop;
        return;
    }

    /* two limbs of padding below both pairs: the window [n, 2n] loses
       up to two limbs against the product's true top when the
       operands' top limbs are partial, and the frame must keep cap */
    n = FLINT_MAX(*ncn, fcn) + 2;
    pn = n - *ncn;
    pf = n - fcn;

    TMP_START;
    ar = TMP_ALLOC(6 * (n + 2) * sizeof(ulong));
    ai = ar + (n + 2);
    br = ai + (n + 2);
    bi = br + (n + 2);
    rr = bi + (n + 2);
    ri = rr + (n + 2);

    flint_mpn_zero(ar, pn);
    flint_mpn_copyi(ar + pn, nc, *ncn);
    flint_mpn_zero(ai, n);
    flint_mpn_copyi(ai + pn, ns, *nsn);
    flint_mpn_zero(br, pf);
    flint_mpn_copyi(br + pf, fc, fcn);
    flint_mpn_zero(bi, n);
    flint_mpn_copyi(bi + pf, fs, fsn);

    flint_mpn_mulhigh_n_complex(rr, &sr, ri, &si, ar, 0, ai, 0, br, 0, bi, 0, n);
    FLINT_ASSERT(sr == 0);

    lr = n + 1;
    while (lr > 1 && rr[lr - 1] == 0)
        lr--;
    li = n + 1;
    while (li > 1 && ri[li - 1] == 0)
        li--;
    if (si)
    {
        /* a vanishing imaginary part rounded below zero */
        ri[0] = 0;
        li = 1;
    }

    /* NUM F = (rr + i ri) B^(nexp - pn + fexp - pf + n); truncate to
       cap limbs keyed to the real part */
    *nexp += fexp - pn - pf + n;
    drop = (lr > cap) ? (lr - cap) : 0;
    *nexp += drop;

    k = FLINT_MAX(1, lr - drop);
    flint_mpn_copyi(nc, rr + drop, k);
    *ncn = k;
    if (li > drop)
    {
        flint_mpn_copyi(ns, ri + drop, li - drop);
        *nsn = li - drop;
    }
    else
    {
        ns[0] = 0;
        *nsn = 1;
    }
    TMP_END;
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
    nn_ptr nc, ns, den, den2, q, ycos;
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
        + (wn + 4)
        + (wn + 1)) * sizeof(ulong));
    ns = nc + (cap + 1);
    den = ns + (cap + 1);
    den2 = den + (cap + 1);
    q = den2 + (cap + 1);
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
            fexp, cap);

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
