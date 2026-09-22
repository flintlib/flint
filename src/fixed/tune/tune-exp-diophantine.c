/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Tune fixed_exp_diophantine: for a size n, sweep the number of
   primes and the weight budget (as a factor of the precision
   FLINT_BITS n), timing the evaluation with a warm cache and
   reporting the residual depth r the reduction reaches; then
   compare the precomputation cost of the diophantine (multi-prime) scheme (one
   logarithm per prime at n limbs) against that of the bitwise
   reduction table at its default depth, and the number of
   evaluations after which the (faster per call) bitwise scheme has
   amortized its (larger) precomputation.

   Usage: tune-exp-diophantine n [nx] */

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "profiler.h"
#include "flint.h"
#include "mpn_extras.h"
#include "arb.h"
#include "fixed.h"

static double
wall_us(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1e6 + tv.tv_usec;
}

/* precomputation of the diophantine (multi-prime) scheme, from scratch: the cache
   build (arb_log_primes_vec_bsplit, floored) */
static double
precomp_diophantine(slong num, slong n)
{
    double t0, t1;
    _fixed_log_primes_clear();
    t0 = wall_us();
    _fixed_log_primes_ensure(num, n + 1);
    t1 = wall_us();
    return t1 - t0;
}

static double
precomp_bitwise(slong n)
{
    double t0, t1;
    nn_ptr y = flint_malloc((n + 1) * sizeof(ulong));
    nn_ptr x = flint_calloc(n, sizeof(ulong));
    x[n - 1] = UWORD(1) << (FLINT_BITS - 2);
    _fixed_exp_logs_clear();
    t0 = wall_us();
    fixed_exp_bitwise_rs(y, x, n, 0);
    t1 = wall_us();
    flint_free(x);
    flint_free(y);
    return t1 - t0;
}

static double
time_call(void (*f)(nn_ptr, nn_srcptr, slong, slong, double),
    nn_ptr y, nn_srcptr xs, slong n, slong nx, slong num, double w)
{
    double best = 1e300;
    int pass;

    for (pass = 0; pass < 2; pass++)
    {
        slong reps = 1;
        timeit_t tm;

        while (1)
        {
            slong it, k;
            timeit_start(tm);
            for (it = 0; it < reps; it++)
                for (k = 0; k < nx; k++)
                    f(y, xs + k * n, n, num, w);
            timeit_stop(tm);
            if (tm->wall >= 40)
                break;
            reps *= 2;
        }
        best = FLINT_MIN(best, tm->wall * 1e3 / (reps * nx));
    }
    return best;
}

static void
f_bitwise(nn_ptr y, nn_srcptr x, slong n, slong num, double w)
{
    (void) num; (void) w;
    fixed_exp_bitwise_rs(y, x, n, 0);
}

/* the residual depth: rerun the reduction and read off the leading
   zeros of the exact residual */
static double
residual_bits(nn_srcptr x, slong n, slong num, double w)
{
    const fixed_rel_struct * tab = fixed_rel_table(0, num);
    slong wr, j;
    slong * rel;
    nn_ptr base, t, tmp;
    double r, eps_min;

    if (FLINT_BITS * n <= 10000) wr = 256 / FLINT_BITS;
    else if (FLINT_BITS * n <= 100000) wr = 512 / FLINT_BITS;
    else wr = 768 / FLINT_BITS;
    wr = FLINT_MAX(wr, (slong) ((-log2(tab->epsilon_min) + 80) / FLINT_BITS) + 1);
    _fixed_log_primes_ensure(num, FLINT_MAX(wr, n + 1));

    rel = flint_malloc(num * sizeof(slong));
    base = flint_malloc(wr * sizeof(ulong));
    t = flint_malloc(2 * (n + 2) * sizeof(ulong));
    tmp = t + (n + 2);
    {
        slong pad = FLINT_MAX(wr - n, 0), i;
        for (i = 0; i < pad; i++) base[i] = 0;
        flint_mpn_copyi(base + pad, x + FLINT_MAX(n - wr, 0), wr - pad);
    }
    eps_min = ldexp(1.0, -(int) FLINT_MIN(FLINT_BITS * n + 32, 2000));
    eps_min = FLINT_MAX(eps_min, ldexp(1.0, -(int) (FLINT_BITS * wr - 48)));
    _fixed_log_reduce(rel, tab, base, wr, w, eps_min,
        _fixed_log_primes_entry(0, wr), _fixed_log_primes_n);

    /* t = x B - sum rel_j log p_j at n + 1 fraction limbs */
    tmp[0] = 0;
    flint_mpn_copyi(tmp + 1, x, n);
    tmp[n + 1] = 0;
    flint_mpn_copyi(t, tmp, n + 2);
    for (j = 0; j < num; j++)
    {
        if (rel[j] > 0)
            mpn_submul_1(t, _fixed_log_primes_entry(j, n + 1), n + 2, rel[j]);
        else if (rel[j] < 0)
            mpn_addmul_1(t, _fixed_log_primes_entry(j, n + 1), n + 2, -rel[j]);
    }
    if (t[n + 1] >> (FLINT_BITS - 1))
        r = -1.0;   /* negative: should not happen */
    else
    {
        slong i;
        for (i = n; i >= 0 && t[i] == 0; i--)
            ;
        r = (i < 0) ? (double) (FLINT_BITS * (n + 1))
            : (double) (FLINT_BITS * (n - i))
              + (double) (FLINT_BITS - (slong) FLINT_BIT_COUNT(t[i]));
    }

    flint_free(rel); flint_free(base); flint_free(t);
    return r;
}

int
main(int argc, char * argv[])
{
    static const slong nums[] = { 4, 6, 8, 10, 13, 16, 20, 24, 32 };
    static const double wfs[] = { 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 32.0 };
    slong n = (argc > 1) ? atol(argv[1]) : 256;
    slong nx = (argc > 2) ? atol(argv[2]) : FLINT_MAX(2, FLINT_MIN(32, 20000 / n));
    slong nnum = sizeof(nums) / sizeof(nums[0]);
    slong nwf = sizeof(wfs) / sizeof(wfs[0]);
    flint_rand_t state;
    nn_ptr xs, y;
    slong k, iw, in;
    double tb, pre_b, pre_m[16], tm_best = 1e300;
    slong best_num = 0;
    double best_wf = 0;

    flint_rand_init(state);
    xs = flint_malloc(nx * n * sizeof(ulong));
    y = flint_malloc((n + 1) * sizeof(ulong));
    for (k = 0; k < nx; k++)
        flint_mpn_urandomb(xs + k * n, state, FLINT_BITS * n);

    /* bitwise reference: precomputation and per-call time */
    pre_b = precomp_bitwise(n);
    pre_b = FLINT_MIN(pre_b, precomp_bitwise(n));
    tb = time_call(f_bitwise, y, xs, n, nx, 0, 0);

    flint_printf("n = %wd (%wd bits), %wd inputs\n", n, FLINT_BITS * n, nx);
    flint_printf("bitwise (r = %d): precomputation %.0f us, per call %.1f us\n\n",
        fixed_exp_bitwise_rs_default_r(n), pre_b, tb);

    flint_printf("per-call time (us) and residual depth r, by primes x weight factor\n");
    flint_printf("primes  precomp(us)");
    for (iw = 0; iw < nwf; iw++)
        flint_printf("      w=%5.2f", wfs[iw]);
    flint_printf("\n");

    for (in = 0; in < nnum; in++)
    {
        slong num = nums[in];

        fixed_rel_table(0, num);   /* table generation outside the timing */
        pre_m[in] = precomp_diophantine(num, n);
        pre_m[in] = FLINT_MIN(pre_m[in], precomp_diophantine(num, n));

        flint_printf("%4wd  %10.0f  ", num, pre_m[in]);
        for (iw = 0; iw < nwf; iw++)
        {
            double w = wfs[iw] * FLINT_BITS * n;
            double t = time_call(_fixed_exp_diophantine_tune, y, xs, n, nx, num, w);
            double r = 0;
            for (k = 0; k < nx; k++)
                r += residual_bits(xs + k * n, n, num, w);
            r /= nx;
            flint_printf(" %7.1f/%3.0f", t, r);
            if (t < tm_best)
            {
                tm_best = t;
                best_num = num;
                best_wf = wfs[iw];
            }
        }
        flint_printf("\n");
    }

    flint_printf("\nbest: %wd primes, w = %.2f x prec: %.1f us per call (%.2fx bitwise)\n",
        best_num, best_wf, tm_best, tm_best / tb);
    for (in = 0; in < nnum; in++)
        if (nums[in] == best_num)
            flint_printf("diophantine precomputation %.0f us vs bitwise %.0f us: "
                "bitwise amortizes its table after %.1f evaluations\n",
                pre_m[in], pre_b, (pre_b - pre_m[in]) / (tm_best - tb));

    flint_free(xs);
    flint_free(y);
    flint_rand_clear(state);
    flint_cleanup();
    return 0;
}
