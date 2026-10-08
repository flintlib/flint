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
#define EXP_REDUCED_MAX_ERR 96

/* _mp_real_exp_reduced: exp(t) of a reduced argument t < 2^-r,
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
   mp_real ball arithmetic (see _mp_real_exp_reduced_ball), whose rigorous
   bound is checked against EXP_REDUCED_MAX_ERR on export (it
   comes to 1 ulp for the one-step cascade and the full burst, and
   the series kernel's own, about 17 ulps, for the sinh path).  The
   series remainder of the cascade runs one guard limb deeper than
   the output: multiplied into the product of the slice factors,
   whose top limb can hold anything from 1 to 64 bits, its radius is
   requantized at a coarser limb anchor, which can triple it per
   cascade level (at most seven levels), and without the guard limb
   that compounded past 96 ulps at four levels (127 ulps at 6292
   limbs, r = 185); with it the remainder's radius stays below one
   ulp of the output at every depth.  An earlier driver in hand-written mpn arithmetic
   (windowed middle products, limb-aligned frames and an informal
   error accounting) measured the same speed at every size and was
   dropped. */

/* Use the sinh series (half the terms) plus a squaring and a square
   root once the direct exp series gets long enough.  Measured
   crossovers on x86-64: about 128 terms (n ~ 2r) at every r below
   384 now that the series tapers at bit granularity (series_rs.c);
   against the earlier untapered series for r < 64 it was 45.  The margins near the crossovers are within a few
   percent, so the exact placement is not critical. */
/* sinh from this many series terms (64 wn / r): the crossover rises
   with r (tune-exp-reduced on the target machine: (64, 128] up to
   r = 256, (171, 256] for r = 384 .. 2048, (256, 512] at 4096) */
/* (at r < 64 the crossover was 45 terms against the untapered series;
   the tapered rectangular splitting of series_rs.c holds to about 128
   there too) */
#define EXP_USE_SINH(wn, r) \
    (FLINT_BITS * (wn) >= (((r) >= 3072) ? 300 : ((r) >= 384) ? 200 \
        : 128) * (slong) (r))

/* When the direct series would need many terms, evaluate sinh(t)
   instead -- the odd series has half the terms -- and reconstruct
   exp(t) = sinh(t) + sqrt(1 + sinh(t)^2), a ball squaring, sum and
   square root on the imported sinh (15 ulp, as returned),
   which comes to about 17 ulps.  One bit-burst step doubles the
   convergence rate: for t < 2^-r with r a whole number of limbs,
   split t = x1 + x2 where x1 is the leading r-bit chunk (bits
   [r, 2r)) and x2 < 2^-2r the remainder, and compute
   exp(t) = exp(x1) exp(x2), exp(x1) = (Q B^QE + T) / (Q B^QE) from
   the mpn binary splitting _mp_real_exp_sum_bs_powtab -- the chunk is
   an r-bit rational, so the splitting tree stays small -- and
   exp(x2) from this function at the doubled rate (possibly bursting
   again). */
/* terms threshold for one burst step; the slice mechanics work for
   any r >= 32 since the boundaries sit on limbs.  Measured
   crossovers (tune-exp-reduced on the target machine): (256, 512]
   at r = 64 and 128, (341, 683] .. (512, 1024] for r = 192 .. 1024,
   (683, 1365] from r = 1536. */
#define MP_REAL_EXP_BURST_TERMS(r) \
    (((r) >= 1536) ? 1000 : ((r) > 128) ? 600 : 512)
/* THREADS.  A burst step's slice runs on as many threads as its
   splitting tree can use, while the series it replaces is one
   indivisible piece (its products use the FFT's threads, at best);
   with t threads the step therefore pays from about 1/t of the
   serial crossover, and the threshold is divided by the thread count
   (capped at MP_REAL_BURST_THREADS_MAX).  The choice thus depends on
   the thread count, as in arb_exp_arf_bb, and so does the exact
   rounding of the result (its accuracy does not). */
#define MP_REAL_BURST_THREADS_MAX 16
#define MP_REAL_BURST_THREADS() \
    FLINT_MIN(flint_get_num_threads(), MP_REAL_BURST_THREADS_MAX)
#define EXP_USE_BURST(wn, r) \
    (FLINT_BITS * (wn) * MP_REAL_BURST_THREADS() \
        >= MP_REAL_EXP_BURST_TERMS(r) * (slong) (r))

/* Full bit-burst from this many terms.  The development VM's value
   (4096) sat below the crossover measured on the target machine at
   every r: (5461, 8192] for r >= 64, (8192, 16384] at r = 32,
   (16384, 32768] at r = 16 -- the one-step cascade, which finishes
   by rectangular splitting, beats the full burst for longer than
   the earlier tuning found.  (The old rationale: one burst step
   from ~512 terms, later than under the per-level divisions whose
   small level-0 divisor was cheaper than the generic rational
   accumulation for one slice.) */
#define MP_REAL_EXP_FULLBURST_TERMS(r) \
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

/* the tuned automatic choice */
static int
_exp_reduced_alg(slong wn, flint_bitcnt_t r)
{
    if (FLINT_BITS * (ulong) wn >= (ulong) MP_REAL_EXP_FULLBURST_TERMS(r) * r)
        return 4;
    else if (r < 32 || EXP_USE_BURST(wn, r))
        return 3;
    else if (EXP_USE_SINH(wn, r))
        return 2;
    else
        return 1;
}

/* ==== the algorithms in mp_real arithmetic ================================

   exp(t) as a ball, with the exponent and error bookkeeping left to
   mp_real: only the two series kernels are imported (_mp_real_exp_rs and
   _mp_real_sinh_rs, with their documented one-sided errors), the sinh
   reconstruction sinh + sqrt(1 + sinh^2) is a ball squaring, sum and
   square root, and in the burst each slice's exact splitting output
   T / (Q B^QEk) becomes the factor Q B^QEk + T by one ball addition
   (which reads only the top cap + 2 limbs of either operand),
   NUM = prod factor and DEN = prod Q are product trees of ball
   products at cap limbs (their operands truncated to the precision
   that contributes), the series remainder of the one-step variant is
   this function again at the rate the cascade ends at, and the finish
   is one ball division.  The bound of the result is rigorous, and
   checked against the documented budget on export (it lands far
   below).

   THREADS.  The slices and the series remainder are independent, so
   they run as tasks on the available threads (_mp_real_parallel_tasks,
   costliest first), the cascade of one-step reductions being
   flattened into the list of slices beforehand so that its levels
   are tasks too rather than a chain of nested calls; the product
   trees fork their halves, and NUM and DEN run side by side.  Above
   MP_REAL_EXP_BURST_SERIAL_LIMBS the slices run one at a time (memory),
   the threads then going to the binary splitting inside each slice
   (exp_sum_bs.c), which forks its subtrees and merges.  This is the
   pair of strategies of arb_exp_arf_bb: parallel exponentials, and
   a serial loop over parallel binary splittings at the largest
   precisions.  The results do not depend on the thread count. */
/* the burst's factors as parallel tasks: tasks 0 .. nb - 2 are the
   slices, in order of decreasing cost (the leading slice has the most
   terms), f[k] = Q B^QE + T of slice k with fq[k] = Q B^QE its
   denominator; a slice of zero bits leaves f[k] = 1.  With series_alg
   nonzero, the last task is exp of the residual below limb depth
   L[nb - 1] by that series algorithm, with fq = 0 meaning no
   denominator.  The series task comes first in the order, being the
   one indivisible piece. */
typedef struct
{
    nn_srcptr t;
    slong wn, cap;
    const slong * L;
    slong nb;
    int series_alg;
    mp_real_struct * f, * fq;
}
exp_burst_struct;

static void _mp_real_exp_reduced_ball(mp_real_t res, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg);

static void
_exp_burst_task(slong i, void * arg)
{
    exp_burst_struct * S = (exp_burst_struct *) arg;
    slong wn = S->wn, cap = S->cap, k;
    nn_srcptr t = S->t;
    const slong * L = S->L;
    mp_real_struct * f, * fq;
    mp_real_phase_t ph;
    TMP_INIT;

    TMP_START;

    if (S->series_alg && i == 0)
    {
        /* the residual below limb depth L[nb - 1] by the series */
        nn_ptr xres;

        k = S->nb - 1;
        f = S->f + k;
        fq = S->fq + k;
        MP_REAL_PHASE_START(ph, "exp burst: series remainder", wn);
        /* one guard limb deeper than the output (a zero limb below
           t): see the notes at the top */
        xres = TMP_ALLOC((wn + 1) * sizeof(ulong));
        xres[0] = 0;
        flint_mpn_copyi(xres + 1, t, wn);
        flint_mpn_zero(xres + 1 + wn - L[k], L[k]);
        _mp_real_exp_reduced_ball(f, xres, wn + 1,
            (flint_bitcnt_t) (FLINT_BITS * L[k]), S->series_alg);
        mp_real_zero(fq);
        MP_REAL_PHASE_END(ph);
    }
    else
    {
        slong D, xn, N, tn, qn, QEk;
        nn_srcptr u;
        nn_ptr T, Q;

        k = S->series_alg ? i - 1 : i;
        f = S->f + k;
        fq = S->fq + k;
        D = L[k + 1];
        u = t + (wn - D);
        xn = D - (k ? L[k] : 0);

        while (xn > 0 && u[xn - 1] == 0)
            xn--;
        if (xn == 0)
        {
            mp_real_set_ui(f, 1);
            mp_real_zero(fq);
            TMP_END;
            return;
        }
        while (xn > 1 && u[0] == 0)
        {
            u++;
            xn--;
            D--;
        }

        N = _mp_real_exp_bs_num_terms((flint_bitcnt_t)
            (FLINT_BITS * D - mpn_bits(u, xn)),
            FLINT_BITS * wn + 64);

        {
            slong t_alloc = N * (D + 2) + 4;
            slong q_alloc2 = (N * FLINT_BIT_COUNT((ulong) N + 1))
                / FLINT_BITS + 3;

            T = TMP_ALLOC((t_alloc + q_alloc2) * sizeof(ulong));
            Q = T + t_alloc;
        }

        if (mp_real_get_verbose() >= 2 && wn >= MP_REAL_VERBOSE_MIN_LIMBS)
            _mp_real_log("exp burst: slice %wd of %wd: %wd limbs at depth %wd, %wd terms",
                k, S->nb - 1, xn, D, N);
        MP_REAL_PHASE_START(ph, "exp burst: slice tree", wn);
        _mp_real_exp_sum_bs_powtab(T, &tn, Q, &qn, &QEk, u, xn, D, N);
        MP_REAL_PHASE_END(ph);
        MP_REAL_PHASE_START(ph, "exp burst: slice factor", wn);

        /* the factor Q B^QEk + T over the denominator Q B^QEk */
        _mp_real_set_mpn_2exp(fq, Q, qn, FLINT_BITS * QEk);
        _mp_real_set_mpn_2exp(f, T, tn, 0);
        mp_real_add(f, f, fq, cap);
        MP_REAL_PHASE_END(ph);
    }

    TMP_END;
}

/* NUM = prod nums and DEN = prod dens as two jobs */
typedef struct
{
    mp_real_struct * num, * den, * nums, * dens;
    slong nn, nd, n;
}
exp_burst_prod_struct;

static void
_exp_burst_prod_num(void * arg)
{
    exp_burst_prod_struct * P = (exp_burst_prod_struct *) arg;
    _mp_real_vec_prod(P->num, P->nums, P->nn, P->n);
}

static void
_exp_burst_prod_den(void * arg)
{
    exp_burst_prod_struct * P = (exp_burst_prod_struct *) arg;
    _mp_real_vec_prod(P->den, P->dens, P->nd, P->n);
}

/* res *= x at n limbs, x destroyed (cleared); a unit res takes x by a
   swap */
static void
_burst_fold(mp_real_t res, mp_real_t x, slong n)
{
    if (res->size == 1 && res->d[0] == 1 && res->exp == 1 && res->err == 0
        && !res->negative)
        mp_real_swap(res, x);
    else
        mp_real_mul(res, res, x, n);
    mp_real_clear(x);
    mp_real_init(x);
}

/* task i in lane w: its factor, folded at once into the lane's running
   products numw[w], denw[w], so that at most two factors per thread
   are ever alive (not one per slice); the lanes and the order within
   each are fixed in advance (_mp_real_parallel_lanes), so the result
   does not depend on timing */
typedef struct
{
    exp_burst_struct * S;
    mp_real_struct * numw, * denw;
}
exp_burst_acc_struct;

static void
_exp_burst_task_acc(slong i, slong w, void * arg)
{
    exp_burst_acc_struct * A = (exp_burst_acc_struct *) arg;
    exp_burst_struct * S = A->S;
    slong j = S->series_alg ? (i ? i - 1 : S->nb - 1) : i;
    mp_real_phase_t ph;

    _exp_burst_task(i, S);
    MP_REAL_PHASE_START(ph, "exp burst: accumulate", S->wn);
    if (!(S->f[j].size == 1 && S->f[j].d[0] == 1 && S->f[j].exp == 1
        && S->f[j].err == 0))
        _burst_fold(A->numw + w, S->f + j, S->cap);
    if (S->fq[j].size != 0)
        _burst_fold(A->denw + w, S->fq + j, S->cap);
    MP_REAL_PHASE_END(ph);
}

static void
_mp_real_exp_reduced_ball(mp_real_t res, nn_srcptr t, slong wn, flint_bitcnt_t r,
    int alg)
{
    slong cap = wn + EXP_BURST_GUARD;

    if (alg == 0)
        alg = _exp_reduced_alg(wn, r);

    if (alg == 1)
    {
        nn_ptr y;
        ulong e;
        TMP_INIT;
        TMP_START;
        y = TMP_ALLOC((wn + 1) * sizeof(ulong));
        _mp_real_exp_rs(y, &e, t, wn);
        _mp_real_set_mpn_2exp(res, y, wn + 1, -FLINT_BITS * wn);
        _mp_real_add_error_ulps_at(res, e, -wn);
        TMP_END;
    }
    else if (alg == 2)
    {
        /* exp = sinh + sqrt(1 + sinh^2) */
        mp_real_t s, u;
        nn_ptr y;
        ulong e;
        TMP_INIT;
        TMP_START;
        y = TMP_ALLOC((wn + 1) * sizeof(ulong));
        mp_real_init(s);
        mp_real_init(u);
        _mp_real_sinh_rs(y, &e, t, wn);
        _mp_real_set_mpn_2exp(s, y, wn + 1, -FLINT_BITS * wn);
        _mp_real_add_error_ulps_at(s, e, -wn);
        mp_real_mul(u, s, s, cap);
        mp_real_set_ui(res, 1);
        mp_real_add(u, u, res, cap);
        mp_real_sqrt(u, u, cap);
        mp_real_add(res, u, s, cap);
        mp_real_clear(s);
        mp_real_clear(u);
        TMP_END;
    }
    else
    {
        slong L[FLINT_BITS + 2];
        slong nb = 0, k, nn, nd, ntasks, nw, maxw;
        int series_alg = 0, mode = alg, serial;
        mp_real_struct * f, * fq, * nums, * dens, * numw, * denw;
        mp_real_t num, den;
        exp_burst_struct S;
        exp_burst_acc_struct A;
        exp_burst_prod_struct P;
        mp_real_phase_t ph;

        /* The ladder of limb depths: slice k covers the limbs between
           depths L[k] and L[k + 1], the depths doubling.  A single step
           (alg 3) is followed by the tuned choice for the residual
           below L[1]; the cascade of steps that choice produces is
           flattened here into the same ladder, ending either at wn
           (the full burst) or with the residual below L[nb - 1] left
           to a series algorithm -- so that all the slices and the
           series remainder can run as tasks side by side. */
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
                mode = _exp_reduced_alg(wn,
                    (flint_bitcnt_t) (FLINT_BITS * L[nb - 1]));
                if (mode <= 2)
                {
                    series_alg = mode;
                    break;
                }
            }
            L[nb] = FLINT_MIN(2 * L[nb - 1], wn);
            nb++;
        }
        if (nb == 1)
        {
            L[1] = wn;
            nb = 2;
        }
        ntasks = (nb - 1) + (series_alg != 0);

        /* The factors, as tasks in order of decreasing cost (the
           series remainder first, being indivisible, then the slices,
           the leading one having the most terms), with cost estimates
           in full products: a slice's tree has about log2 N levels,
           the series about 2 sqrt(N) products.  The tasks are dealt
           to one lane per thread by estimated cost, and each lane
           folds its factors into its own running numerator and
           denominator as they complete, so that the live
           full-precision numbers are two per lane plus the ones being
           computed, independently of the number of slices; the lanes'
           products are then combined by product trees.  In the serial mode -- above
           MP_REAL_EXP_BURST_SERIAL_LIMBS, where a slice's splitting
           integers are large enough for memory to matter, and below
           MP_REAL_PAR_MIN_LIMBS, where a thread is not worth waking --
           a single worker runs the slices one at a time, the threads
           going to the binary splitting and the products inside. */
        serial = (wn > MP_REAL_EXP_BURST_SERIAL_LIMBS || wn < MP_REAL_PAR_MIN_LIMBS);
        maxw = serial ? 1 : ntasks;
        nw = _mp_real_parallel_lanes_count(ntasks, maxw);

        f = flint_malloc((2 * nb + 4 * nw) * sizeof(mp_real_struct));
        fq = f + nb;
        numw = fq + nb;
        denw = numw + nw;
        nums = denw + nw;
        dens = nums + nw;
        for (k = 0; k < 2 * nb + 4 * nw; k++)
            mp_real_init(f + k);
        for (k = 0; k < nw; k++)
        {
            mp_real_set_ui(numw + k, 1);
            mp_real_set_ui(denw + k, 1);
        }
        mp_real_init(num);
        mp_real_init(den);

        S.t = t;
        S.wn = wn;
        S.cap = cap;
        S.L = L;
        S.nb = nb;
        S.series_alg = series_alg;
        S.f = f;
        S.fq = fq;
        A.S = &S;
        A.numw = numw;
        A.denw = denw;
        if (mp_real_get_verbose() >= 2 && wn >= MP_REAL_VERBOSE_MIN_LIMBS)
            _mp_real_log("exp burst: %wd limbs, r = %wd, %wd slices%s, %wd workers",
                wn, (slong) r, nb - 1, series_alg ? " and a series remainder" : "",
                nw);

        {
            double cost[FLINT_BITS + 2];
            slong j = 0;

            if (series_alg)
                cost[j++] = 2.0 * sqrt((double) wn / (double) L[nb - 1]) + 4.0;
            for (k = 0; k < nb - 1; k++)
                cost[j++] = log((double) wn / (double) L[k] + 2.0) * 1.5 + 2.0;

            _mp_real_parallel_lanes(_exp_burst_task_acc, &A, ntasks, nw, cost);
        }

        /* NUM = prod numw, DEN = prod denw over the lanes (skipping the unit ones), the
           two trees on two threads */
        nn = nd = 0;
        for (k = 0; k < nw; k++)
        {
            if (!(numw[k].size == 1 && numw[k].d[0] == 1 && numw[k].exp == 1
                && numw[k].err == 0))
                mp_real_swap(nums + nn++, numw + k);
            if (!(denw[k].size == 1 && denw[k].d[0] == 1 && denw[k].exp == 1
                && denw[k].err == 0))
                mp_real_swap(dens + nd++, denw + k);
        }
        P.num = num; P.den = den; P.nums = nums; P.dens = dens;
        P.nn = nn; P.nd = nd; P.n = cap;
        MP_REAL_PHASE_START(ph, "exp burst: product trees", wn);
        if (nn >= 2 && nd >= 2 && flint_get_num_threads() >= 2)
            _mp_real_parallel_pair(_exp_burst_prod_num, &P,
                _exp_burst_prod_den, &P);
        else
        {
            _exp_burst_prod_num(&P);
            _exp_burst_prod_den(&P);
        }
        MP_REAL_PHASE_END(ph);

        MP_REAL_PHASE_START(ph, "exp burst: final division", wn);
        mp_real_div(res, num, den, wn + 2);
        MP_REAL_PHASE_END(ph);

        for (k = 0; k < 2 * nb + 4 * nw; k++)
            mp_real_clear(f + k);
        mp_real_clear(num);
        mp_real_clear(den);
        flint_free(f);
    }
}

/* the mp_real algorithms exported in the fixed-point format */
static void
_mp_real_exp_reduced_export(nn_ptr y, ulong * err, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg)
{
    mp_real_t v, one;
    ulong bound;

    /* the plain series is already in the export format: skip the ball
       conversions (their fixed cost is comparable to the whole series
       at small sizes) */
    if (alg == 1 || (alg == 0 && _exp_reduced_alg(wn, r) == 1))
    {
        _mp_real_exp_rs(y, err, t, wn);
        return;
    }

    mp_real_init(v);
    mp_real_init(one);
    _mp_real_exp_reduced_ball(v, t, wn, r, alg);
    /* exp(t) - 1 in [0, 1) */
    mp_real_set_ui(one, 1);
    mp_real_sub(v, v, one, wn + 2);
    _mp_real_get_fixed(y, &bound, v, wn);
    if (!(bound < EXP_REDUCED_MAX_ERR))
        flint_throw(FLINT_ERROR, "_mp_real_exp_reduced (mp_real): error bound %wu ulps\n", bound);
    y[wn] = 1;
    if (err != NULL)
        *err = bound;
    mp_real_clear(v);
    mp_real_clear(one);
}

void
_mp_real_exp_reduced(nn_ptr y, ulong * err, nn_srcptr t, slong n,
    flint_bitcnt_t r, int alg)
{
    FLINT_ASSERT(r >= 16);
    FLINT_ASSERT(alg == 0 || alg >= 3 || r >= 32);
    _mp_real_exp_reduced_export(y, err, t, n, r, alg);
}
