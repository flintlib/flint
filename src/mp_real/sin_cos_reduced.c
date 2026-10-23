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
   Instead of combining per level, the complex numerator
   NUM = prod (Q_k B^{QE_k} - A_k + i B_k) and the real denominator
   DEN = prod Q_k are product trees of balls at wn + 3 limbs (one
   transform-sharing complex product per node, mp_real_mul_complex),
   and the finish is TWO ball divisions (sine and cosine) against the
   single accumulated denominator -- where arb_sin_cos_arf_bb spends one full-precision
   square root PER SLICE (cos_k from sin_k) plus a full-precision
   complex product tree.  g comes from the cosine by subtraction:
   the ~2r-bit cancellation is harmless because g is only needed to
   absolute 2^(-64 wn), like the direct series path.  The bounds are
   rigorous, checked against SIN_COS_REDUCED_MAX_ERR = 96 ulps
   on export: 1 ulp for the burst variants, and the series kernels'
   own for the sine + sqrt path.  The series remainder of the
   one-step variant runs one guard limb deeper than the output, as in
   _mp_real_exp_reduced: multiplied into the product of the slice
   factors, its radius is requantized at a coarser limb anchor
   (up to three times per level), which without the guard limb
   compounded to 70 ulps (5901 limbs, r = 31).  An earlier driver in
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
   terms per r (8.4 10^6 bits at r = 16), nor at 10^8 bits (63 against
   55 seconds for the cascade on the development VM, the deepest
   slices' trees and factors costing more than the series they stand
   in for), so the cascade is always used; the constant remains for
   experiments */
#ifndef MP_REAL_TRIG_FULLBURST_TERMS
#if FLINT_BITS == 64
#define MP_REAL_TRIG_FULLBURST_TERMS (WORD(1) << 40)
#else
/* (out of reach as well; a shift by 40 would be undefined) */
#define MP_REAL_TRIG_FULLBURST_TERMS (WORD(1) << 26)
#endif
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
/* with threads, from 1/t of the serial crossover, as in exp_reduced.c */
#define MP_REAL_BURST_THREADS_MAX 16
#define MP_REAL_BURST_THREADS() \
    FLINT_MIN(flint_get_num_threads(), MP_REAL_BURST_THREADS_MAX)
#define TRIG_USE_BURST(wn, r) \
    (FLINT_BITS * (wn) * MP_REAL_BURST_THREADS() \
        >= (((r) < 64) ? MP_REAL_TRIG_BURST_TERMS_SMALL_R \
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
   is this function again at the rate the cascade ends at, and the
   finish is two ball divisions.  The bounds are rigorous and checked
   against the documented budget on export.  The slices and the
   series remainder run as parallel tasks, the product trees fork,
   and above MP_REAL_TRIG_BURST_SERIAL_LIMBS the slices go serial with the
   threads inside their binary splitting (sin_cos_sum_bs.c), exactly
   as in _mp_real_exp_reduced_ball. */
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

/* the burst's factors as parallel tasks: tasks 0 .. nb - 2 are the
   slices, in order of decreasing cost (the leading slice has the most
   terms), fc[k] + i fs[k] the factor of slice k with fq[k] = Q B^QE its
   denominator; a slice of zero bits leaves fc[k] = 1, fs[k] = 0.  With
   series_alg nonzero, the last task is cos + i sin of the residual
   below limb depth L[nb - 1] by that series algorithm, with fq = 0
   meaning no denominator.  The series task comes first in the order,
   being the one indivisible piece. */
typedef struct
{
    nn_srcptr t;
    slong wn, cap;
    const slong * L;
    slong nb;
    int series_alg;
    mp_real_struct * fc, * fs, * fq;
}
trig_burst_struct;

static void _mp_real_sin_cos_reduced_ball(mp_real_t rs, mp_real_t rc,
    nn_srcptr t, slong wn, flint_bitcnt_t r, int alg);

static void
_trig_burst_task(slong i, void * arg)
{
    trig_burst_struct * S = (trig_burst_struct *) arg;
    slong wn = S->wn, cap = S->cap, k;
    nn_srcptr t = S->t;
    const slong * L = S->L;
    mp_real_struct * fc, * fs, * fq;
    mp_real_phase_t ph;
    TMP_INIT;

    TMP_START;

    if (S->series_alg && i == 0)
    {
        nn_ptr xres;

        k = S->nb - 1;
        fc = S->fc + k;
        fs = S->fs + k;
        fq = S->fq + k;
        MP_REAL_PHASE_START(ph, "trig burst: series remainder", wn);
        /* one guard limb deeper than the output (a zero limb below
           t): see the notes at the top */
        xres = TMP_ALLOC((wn + 1) * sizeof(ulong));
        xres[0] = 0;
        flint_mpn_copyi(xres + 1, t, wn);
        flint_mpn_zero(xres + 1 + wn - L[k], L[k]);
        _mp_real_sin_cos_reduced_ball(fs, fc, xres, wn + 1,
            (flint_bitcnt_t) (FLINT_BITS * L[k]), S->series_alg);
        mp_real_zero(fq);
        MP_REAL_PHASE_END(ph);
    }
    else
    {
        slong D, xn;
        nn_srcptr x;
        slong N, an, bn2, qn, ae, be, QEk;
        nn_ptr A, B, Q;
        int slice_sqrt;
        mp_real_t u, v;

        k = S->series_alg ? i - 1 : i;
        fc = S->fc + k;
        fs = S->fs + k;
        fq = S->fq + k;
        D = L[k + 1];
        x = t + (wn - D);
        xn = D - (k ? L[k] : 0);

        while (xn > 0 && x[xn - 1] == 0)
            xn--;
        if (xn == 0)
        {
            mp_real_set_ui(fc, 1);
            mp_real_zero(fs);
            mp_real_zero(fq);
            TMP_END;
            return;
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
        if (mp_real_get_verbose() >= 2 && wn >= MP_REAL_VERBOSE_MIN_LIMBS)
            _mp_real_log("trig burst: slice %wd of %wd: %wd limbs at depth %wd, %wd terms%s",
                k, S->nb - 1, xn, D, N, slice_sqrt ? ", sine track only" : "");
        MP_REAL_PHASE_START(ph, "trig burst: slice tree", wn);
        _mp_real_sin_cos_sum_bs_powtab(
            slice_sqrt ? NULL : A, &an, &ae, B, &bn2, &be,
            Q, &qn, &QEk, x, xn, D, N, cap + 2);
        MP_REAL_PHASE_END(ph);
        MP_REAL_PHASE_START(ph, "trig burst: slice factor", wn);

        mp_real_init(u);
        mp_real_init(v);

        /* the denominator Q B^QEk of the slice */
        _mp_real_set_mpn_2exp(fq, Q, qn, FLINT_BITS * QEk);

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
        MP_REAL_PHASE_END(ph);

        mp_real_clear(u);
        mp_real_clear(v);
    }

    TMP_END;
}

/* NUM = prod (fc + i fs) and DEN = prod fq as two jobs */
typedef struct
{
    mp_real_struct * nc, * ns, * den, * fc, * fs, * dens;
    slong nn, nd, n;
}
trig_burst_prod_struct;

static void
_trig_burst_prod_num(void * arg)
{
    trig_burst_prod_struct * P = (trig_burst_prod_struct *) arg;
    _mp_real_vec_prod_complex(P->nc, P->ns, P->fc, P->fs, P->nn, P->n);
}

static void
_trig_burst_prod_den(void * arg)
{
    trig_burst_prod_struct * P = (trig_burst_prod_struct *) arg;
    _mp_real_vec_prod(P->den, P->dens, P->nd, P->n);
}

/* task i in lane w, its factor folded at once into the lane's running
   products (see exp_reduced.c) */
typedef struct
{
    trig_burst_struct * S;
    mp_real_struct * ncw, * nsw, * denw;
}
trig_burst_acc_struct;

static int
_is_one(const mp_real_t x)
{
    return x->size == 1 && x->d[0] == 1 && x->exp == 1 && x->err == 0
        && !x->negative;
}

static void
_trig_burst_task_acc(slong i, slong w, void * arg)
{
    trig_burst_acc_struct * A = (trig_burst_acc_struct *) arg;
    trig_burst_struct * S = A->S;
    slong j = S->series_alg ? (i ? i - 1 : S->nb - 1) : i;
    mp_real_struct * fc = S->fc + j, * fs = S->fs + j, * fq = S->fq + j;
    mp_real_phase_t ph;

    _trig_burst_task(i, S);
    MP_REAL_PHASE_START(ph, "trig burst: accumulate", S->wn);
    if (!(_is_one(fc) && fs->size == 0 && fs->err == 0))
    {
        if (_is_one(A->ncw + w) && A->nsw[w].size == 0 && A->nsw[w].err == 0)
        {
            mp_real_swap(A->ncw + w, fc);
            mp_real_swap(A->nsw + w, fs);
        }
        else
            mp_real_mul_complex(A->ncw + w, A->nsw + w, A->ncw + w,
                A->nsw + w, fc, fs, S->cap);
    }
    if (fq->size != 0)
    {
        if (_is_one(A->denw + w))
            mp_real_swap(A->denw + w, fq);
        else
            mp_real_mul(A->denw + w, A->denw + w, fq, S->cap);
    }
    mp_real_clear(fc);
    mp_real_init(fc);
    mp_real_clear(fs);
    mp_real_init(fs);
    mp_real_clear(fq);
    mp_real_init(fq);
    MP_REAL_PHASE_END(ph);
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
        mp_real_phase_t ph;
        mp_real_init(u);
        MP_REAL_PHASE_START(ph, "trig series: sine", wn);
        _mp_real_sin_rs(ss, &e, t, wn);
        MP_REAL_PHASE_END(ph);
        _mp_real_set_mpn_2exp(rs, ss, wn + 1, -FLINT_BITS * wn);
        _mp_real_add_error_ulps_at(rs, e, -wn);
        MP_REAL_PHASE_START(ph, "trig series: square", wn);
        mp_real_mul(u, rs, rs, cap);
        MP_REAL_PHASE_END(ph);
        mp_real_set_ui(rc, 1);
        mp_real_sub(u, rc, u, cap);
        MP_REAL_PHASE_START(ph, "trig series: square root", wn);
        mp_real_sqrt(rc, u, cap);
        MP_REAL_PHASE_END(ph);
        mp_real_clear(u);
        TMP_END;
    }
    else
    {
        slong L[FLINT_BITS + 2];
        slong nb = 0, k, nn, nd, ntasks, nw, maxw;
        int series_alg = 0, mode = alg, serial;
        mp_real_struct * fc, * fs, * fq, * vc, * vs, * dens, * ncw, * nsw, * denw;
        mp_real_t nc, ns, den;
        trig_burst_struct S;
        trig_burst_acc_struct A;
        trig_burst_prod_struct P;
        mp_real_phase_t ph;

        /* the ladder of limb depths (tripling), the one-step cascade
           flattened into it as in _mp_real_exp_reduced_ball */
        L[nb++] = FLINT_MAX((slong) r / FLINT_BITS, 1);
        /* Below a limb of reduction the first slice is the top limb
           alone: a slice's splitting integers grow by its frame (64 D
           bits, 128 D for the sine and cosine's x^2) per term while its
           argument decays by r bits per term, and the ladder keeps
           that ratio at 2 (3 for the tripling ladder) except for a
           first slice of 2 resp. 3 limbs at r < 64, where it is
           128 / r resp. 384 / (2 r) -- at r = 32 the leading slice then
           held integers twice (resp. three times) as large as needed,
           and was both the costliest slice and the peak of the memory. */
        if ((slong) r < FLINT_BITS && wn > 1)
            L[nb++] = 1;
        while (L[nb - 1] < wn)
        {
            if (nb >= 2 && mode == 3)
            {
                mode = _sin_cos_reduced_alg(wn,
                    (flint_bitcnt_t) (FLINT_BITS * L[nb - 1]));
                if (mode <= 2)
                {
                    series_alg = mode;
                    break;
                }
            }
            L[nb] = FLINT_MIN(3 * L[nb - 1], wn);
            nb++;
        }
        if (nb == 1)
        {
            L[1] = wn;
            nb = 2;
        }
        ntasks = (nb - 1) + (series_alg != 0);

        /* the lanes, budgets, running products and the serial mode as
           in _mp_real_exp_reduced_ball */
        serial = (wn > MP_REAL_TRIG_BURST_SERIAL_LIMBS || wn < MP_REAL_PAR_MIN_LIMBS);
        maxw = serial ? 1 : ntasks;
        nw = _mp_real_parallel_lanes_count(ntasks, maxw);

        fc = flint_malloc((3 * nb + 6 * nw) * sizeof(mp_real_struct));
        fs = fc + nb;
        fq = fs + nb;
        ncw = fq + nb;
        nsw = ncw + nw;
        denw = nsw + nw;
        vc = denw + nw;
        vs = vc + nw;
        dens = vs + nw;
        for (k = 0; k < 3 * nb + 6 * nw; k++)
            mp_real_init(fc + k);
        for (k = 0; k < nw; k++)
        {
            mp_real_set_ui(ncw + k, 1);
            mp_real_set_ui(denw + k, 1);
        }
        mp_real_init(nc);
        mp_real_init(ns);
        mp_real_init(den);

        S.t = t;
        S.wn = wn;
        S.cap = cap;
        S.L = L;
        S.nb = nb;
        S.series_alg = series_alg;
        S.fc = fc;
        S.fs = fs;
        S.fq = fq;
        A.S = &S;
        A.ncw = ncw;
        A.nsw = nsw;
        A.denw = denw;
        if (mp_real_get_verbose() >= 2 && wn >= MP_REAL_VERBOSE_MIN_LIMBS)
            _mp_real_log("trig burst: %wd limbs, r = %wd, %wd slices%s, %wd workers",
                wn, (slong) r, nb - 1, series_alg ? " and a series remainder" : "",
                nw);

        {
            double cost[FLINT_BITS + 2];
            slong j = 0;

            if (series_alg)
                cost[j++] = 2.0 * sqrt((double) wn / (2.0 * (double) L[nb - 1])) + 4.0;
            for (k = 0; k < nb - 1; k++)
                cost[j++] = log((double) wn / (2.0 * (double) L[k]) + 2.0) * 2.0 + 3.0;

            _mp_real_parallel_lanes(_trig_burst_task_acc, &A, ntasks, nw, cost);
        }

        /* NUM = prod (ncw + i nsw), DEN = prod denw (skipping the unit
           ones), the two trees on two threads */
        nn = nd = 0;
        for (k = 0; k < nw; k++)
        {
            if (!(_is_one(ncw + k) && nsw[k].size == 0 && nsw[k].err == 0))
            {
                mp_real_swap(vc + nn, ncw + k);
                mp_real_swap(vs + nn, nsw + k);
                nn++;
            }
            if (!_is_one(denw + k))
                mp_real_swap(dens + nd++, denw + k);
        }
        P.nc = nc; P.ns = ns; P.den = den; P.fc = vc; P.fs = vs;
        P.dens = dens; P.nn = nn; P.nd = nd; P.n = cap;
        MP_REAL_PHASE_START(ph, "trig burst: product trees", wn);
        if (nn >= 2 && nd >= 2 && flint_get_num_threads() >= 2)
            _mp_real_parallel_pair(_trig_burst_prod_num, &P,
                _trig_burst_prod_den, &P);
        else
        {
            _trig_burst_prod_num(&P);
            _trig_burst_prod_den(&P);
        }
        MP_REAL_PHASE_END(ph);

        MP_REAL_PHASE_START(ph, "trig burst: final divisions", wn);
        mp_real_div(rs, ns, den, wn + 2);
        mp_real_div(rc, nc, den, wn + 2);
        MP_REAL_PHASE_END(ph);

        for (k = 0; k < 3 * nb + 6 * nw; k++)
            mp_real_clear(fc + k);
        mp_real_clear(nc);
        mp_real_clear(ns);
        mp_real_clear(den);
        flint_free(fc);
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
