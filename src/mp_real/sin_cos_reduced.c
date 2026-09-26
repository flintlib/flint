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
#include "mp_real.h"
#include "impl.h"

/* the documented maximum of *err, checked on export */
#define SIN_COS_REDUCED_MAX_ERR 96

/* _mp_real_sin_cos_reduced: sin(t) and g = 1 - cos(t) of a reduced
   argument t < 2^-r, r >= 16 (the series algorithms 1 and 2 require
   r >= 32), into (ysin, wn) and (yg, wn) -- wn fraction limbs each;
   sin(t) < 2^-r and g < 2^(-2r-1), so neither carries a units limb,
   though both buffers must have room for wn + 1 limbs (the top limb
   is scratch).  This is the trigonometric analogue of
   _mp_real_exp_reduced, independent of any particular
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

   The burst mirrors _mp_real_exp_reduced at LIMB granularity: the
   boundary ladder is a list of limb counts (tripling), a slice is a
   pointer + length into t with a whole-limb mask, dead low limbs of
   a slice fold into its frame (u B^-D = (u/B^z) B^-(D-z)), and each
   slice's factor is the exact Gaussian rational

       exp(i x_k) ~ (Q_k B^{QE_k} - A_k + i B_k) / (Q_k B^{QE_k})

   truncated at N_k series terms, all exponents limb counts.
   Instead of combining per level, the loop accumulates the complex
   numerator NUM = prod (Q_k B^{QE_k} - A_k + i B_k) and the real
   denominator DEN = prod Q_k as balls at wn + 3 limbs (one
   transform-sharing complex product per slice, mp_real_mul_complex),
   deepest slice first so every accumulation product is balanced
   against the content gathered so far, and finishes with TWO ball
   divisions (sine and cosine) against the single accumulated
   denominator -- where arb_sin_cos_arf_bb spends one full-precision
   square root PER SLICE (cos_k from sin_k) plus a full-precision
   complex product tree.  g comes from the cosine by subtraction:
   the ~2r-bit cancellation is harmless because g is only needed to
   absolute 2^(-64 wn), like the direct series path.  The bounds are
   rigorous, checked against SIN_COS_REDUCED_MAX_ERR = 96 ulps
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
#define MP_REAL_TRIG_REDUCED_SINSQRT_TERMS(r) \
    (((r) >= 4096) ? 200 : ((r) >= 2048) ? 100 : 56)
/* r < 64: 320 against the untapered series; with the tapered
   rectangular splitting of series_rs.c the sine + sqrt path holds to
   768 terms at r = 32 and the burst step wins from about 1024 */
#ifndef MP_REAL_TRIG_BURST_TERMS_SMALL_R
#define MP_REAL_TRIG_BURST_TERMS_SMALL_R 896
#endif
#ifndef MP_REAL_TRIG_BURST_TERMS
#define MP_REAL_TRIG_BURST_TERMS 1500
#endif
/* the full bit-burst never won on the target machine up to 524288
   terms per r (8.4 10^6 bits at r = 16); the one-step cascade covers
   that range, so the switch sits above it */
#ifndef MP_REAL_TRIG_FULLBURST_TERMS
#define MP_REAL_TRIG_FULLBURST_TERMS 1048576
#endif
/* Per-slice choice inside the burst: from this many series terms
   the slice's 1 - cos track (two of the four heavy tree
   multiplications per merge) costs more than recovering the cosine
   window by a fixed-point square root -- one unbalanced division
   by the slice denominator, a squaring, the square root and one
   short product against Q.  Below it the joint tree wins (the
   square-root path's costs are flat in N). */
#ifndef MP_REAL_TRIG_SLICE_SQRT_TERMS
#define MP_REAL_TRIG_SLICE_SQRT_TERMS 64
#endif
/* with the tapered series (series_rs.c) the joint sine and cosine
   series stay ahead up to about 48 terms at small sizes, but from 48
   limbs on the square root is cheap enough relative to the second
   series that the sine + sqrt path wins from about 24 terms (measured
   at r = 32 .. 512, wn = 16 .. 2048) */
#define TRIG_USE_SINSQRT(wn, r) \
    (FLINT_BITS * (wn) >= MP_REAL_TRIG_REDUCED_SINSQRT_TERMS(r) * (slong) (r) \
     || ((wn) >= 48 && (r) < 2048 && FLINT_BITS * (wn) >= 24 * (slong) (r)))
#define TRIG_USE_BURST(wn, r) \
    (FLINT_BITS * (wn) >= (((r) < 64) ? MP_REAL_TRIG_BURST_TERMS_SMALL_R \
        : MP_REAL_TRIG_BURST_TERMS) * (slong) (r))


/* ==== the algorithms in mp_real arithmetic ================================

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
/* the tuned automatic choice */
static int
_sin_cos_reduced_alg(slong wn, flint_bitcnt_t r)
{
    ulong terms = FLINT_BITS * (ulong) wn;
    if (terms >= (ulong) MP_REAL_TRIG_FULLBURST_TERMS * r)
        return 4;
    else if (r < 32 || TRIG_USE_BURST(wn, r))
        return 3;
    else if (TRIG_USE_SINSQRT(wn, r))
        return 2;
    else
        return 1;
}

static void
_mp_real_sin_cos_reduced_ball(mp_real_t rs, mp_real_t rc, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg)
{
    slong cap = wn + TRIG_BURST_GUARD;

    if (alg == 0)
        alg = _sin_cos_reduced_alg(wn, r);

    if (alg == 1)
    {
        nn_ptr ss, cc;
        ulong e;
        TMP_INIT;
        TMP_START;
        ss = TMP_ALLOC(2 * (wn + 2) * sizeof(ulong));
        cc = ss + (wn + 2);
        _mp_real_sin_cos_rs(ss, cc, &e, t, wn);
        _mp_real_set_mpn_2exp(rs, ss, wn + 1, -FLINT_BITS * wn);
        _mp_real_add_error_ulps_at(rs, e, -wn);
        _mp_real_set_mpn_2exp(rc, cc, wn + 1, -FLINT_BITS * wn);
        _mp_real_add_error_ulps_at(rc, e, -wn);
        TMP_END;
    }
    else if (alg == 2)
    {
        /* cos = sqrt(1 - sin^2) */
        mp_real_t u;
        nn_ptr ss;
        ulong e;
        TMP_INIT;
        TMP_START;
        ss = TMP_ALLOC((wn + 2) * sizeof(ulong));
        mp_real_init(u);
        _mp_real_sin_rs(ss, &e, t, wn);
        _mp_real_set_mpn_2exp(rs, ss, wn + 1, -FLINT_BITS * wn);
        _mp_real_add_error_ulps_at(rs, e, -wn);
        mp_real_mul(u, rs, rs, cap);
        mp_real_set_ui(rc, 1);
        mp_real_sub(u, rc, u, cap);
        mp_real_sqrt(rc, u, cap);
        mp_real_clear(u);
        TMP_END;
    }
    else
    {
        int levels = (alg == 4) ? 0 : 1;
        slong L[FLINT_BITS + 2];
        slong nb = 0, k;
        mp_real_t nc, ns, den, fc, fs, fq, u, v;

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

        mp_real_init(nc);
        mp_real_init(ns);
        mp_real_init(den);
        mp_real_init(fc);
        mp_real_init(fs);
        mp_real_init(fq);
        mp_real_init(u);
        mp_real_init(v);
        mp_real_set_ui(nc, 1);
        mp_real_set_ui(den, 1);

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
                _mp_real_sin_cos_reduced_ball(fs, fc, xres, wn,
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
                    N = _mp_real_exp_bs_num_terms(
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

                slice_sqrt = (N >= MP_REAL_TRIG_SLICE_SQRT_TERMS);
                _mp_real_sin_cos_sum_bs_powtab(
                    slice_sqrt ? NULL : A, &an, &ae, B, &bn2, &be,
                    Q, &qn, &QEk, x, xn, D, N, cap + 2);

                /* the denominator Q B^QEk of the slice */
                _mp_real_set_mpn_2exp(fq, Q, qn, FLINT_BITS * QEk);
                mp_real_mul(den, den, fq, cap);

                /* FS = x (Q B^QEk - B B^be) B^-D */
                _mp_real_set_mpn_2exp(fs, B, bn2, FLINT_BITS * be);
                mp_real_sub(fs, fq, fs, cap);
                _mp_real_set_mpn_2exp(u, x, xn, -FLINT_BITS * D);
                mp_real_mul(fs, fs, u, cap);

                if (!slice_sqrt)
                {
                    /* FC = Q B^QEk - A B^ae */
                    _mp_real_set_mpn_2exp(fc, A, an, FLINT_BITS * ae);
                    mp_real_sub(fc, fq, fc, cap);
                }
                else
                {
                    /* FC = sqrt((Q B^QEk)^2 - FS^2) */
                    mp_real_mul(u, fq, fq, cap);
                    mp_real_mul(v, fs, fs, cap);
                    mp_real_sub(u, u, v, cap);
                    mp_real_sqrt(fc, u, cap);
                }
            }

            /* NUM *= FC + i FS */
            if (nc->size == 1 && nc->d[0] == 1 && nc->err == 0
                && nc->exp == 1 && ns->size == 0 && ns->err == 0)
            {
                mp_real_set(nc, fc);
                mp_real_set(ns, fs);
            }
            else
                mp_real_mul_complex(nc, ns, nc, ns, fc, fs, cap);
            TMP_END;
        }

        mp_real_div(rs, ns, den, wn + 2);
        mp_real_div(rc, nc, den, wn + 2);

        mp_real_clear(nc);
        mp_real_clear(ns);
        mp_real_clear(den);
        mp_real_clear(fc);
        mp_real_clear(fs);
        mp_real_clear(fq);
        mp_real_clear(u);
        mp_real_clear(v);
    }
}

/* the mp_real algorithms exported in the fixed-point format */
static void
_mp_real_sin_cos_reduced_export(nn_ptr ysin, nn_ptr yg, ulong * err,
    nn_srcptr t, slong wn, flint_bitcnt_t r, int alg)
{
    mp_real_t s, c, one;
    ulong bound, bound2;

    /* the plain series: straight to the fixed-point format, skipping
       the ball conversions (their fixed cost is comparable to the whole
       series at small sizes); g = 1 - cos exactly, since cos <= 1 */
    if (alg == 1 || (alg == 0 && _sin_cos_reduced_alg(wn, r) == 1))
    {
        _mp_real_sin_cos_rs(ysin, yg, err, t, wn);
        if (yg[wn] != 0)
            flint_mpn_zero(yg, wn);         /* cos = 1 */
        else
            mpn_neg(yg, yg, wn);
        return;
    }

    mp_real_init(s);
    mp_real_init(c);
    mp_real_init(one);
    _mp_real_sin_cos_reduced_ball(s, c, t, wn, r, alg);
    _mp_real_get_fixed(ysin, &bound, s, wn);
    if (!(bound < SIN_COS_REDUCED_MAX_ERR))
        flint_throw(FLINT_ERROR, "_mp_real_sin_cos_reduced (mp_real): sine error bound %wu ulps\n", bound);
    /* g = 1 - cos in [0, 1) */
    mp_real_set_ui(one, 1);
    mp_real_sub(c, one, c, wn + 2);
    _mp_real_get_fixed(yg, &bound2, c, wn);
    if (!(bound2 < SIN_COS_REDUCED_MAX_ERR))
        flint_throw(FLINT_ERROR, "_mp_real_sin_cos_reduced (mp_real): cosine error bound %wu ulps\n", bound2);
    if (err != NULL)
        *err = FLINT_MAX(bound, bound2);
    mp_real_clear(s);
    mp_real_clear(c);
    mp_real_clear(one);
}

void
_mp_real_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, ulong * err, nn_srcptr t,
    slong n, flint_bitcnt_t r, int alg)
{
    FLINT_ASSERT(r >= 16);
    FLINT_ASSERT(alg == 0 || alg >= 3 || r >= 32);
    _mp_real_sin_cos_reduced_export(ysin, yg, err, t, n, r, alg);
}
