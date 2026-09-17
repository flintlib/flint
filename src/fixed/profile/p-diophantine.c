/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Compares the argument-reduction strategies for exp and sin_cos at
   a range of precisions: the table-free routines (fixed_exp_notab,
   fixed_sin_cos_notab), the bitwise reductions with their tuned
   table depth r (fixed_exp_bitwise_rs, fixed_sin_cos_bitwise_rs),
   and the diophantine reductions with several prime counts
   (_fixed_exp_diophantine_tune, _fixed_sin_cos_diophantine_tune).

   For each precision the table lists, in multiples of the notab
   time at that precision,

       precomp   the first-call cost beyond a warm call: the bitwise
                 routines' r-bit table, or the diophantine routines'
                 log(p) / angle cache (the relation tables are
                 precomputed for the prime counts used here);
       warm      the evaluation time with the caches filled,

   the tuned r of the bitwise routines, and the ratio of the two
   notab times.  Times are minima over a few random arguments.

   Usage: p-diophantine [n1 n2 ...]     (limbs; default 256 .. 16384)
          -w W   weight budget factor for the diophantine methods
                 (times the precision; default 4)
          -r R   repetitions per argument (default 2) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

#define NX 3

static double
wall(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

typedef void (*exp_fn)(nn_ptr, nn_srcptr, slong, slong, double);
typedef void (*sc_fn)(nn_ptr, nn_ptr, nn_srcptr, slong, slong, double);

static void e_notab(nn_ptr y, nn_srcptr x, slong n, slong num, double w)
{ (void) num; (void) w; fixed_exp_notab(y, x, n); }
static void e_bitwise(nn_ptr y, nn_srcptr x, slong n, slong num, double w)
{ (void) num; (void) w; fixed_exp_bitwise_rs(y, x, n, 0); }
static void e_dio(nn_ptr y, nn_srcptr x, slong n, slong num, double w)
{ _fixed_exp_diophantine_tune(y, x, n, num, w); }
static void s_notab(nn_ptr s, nn_ptr c, nn_srcptr x, slong n, slong num, double w)
{ (void) num; (void) w; fixed_sin_cos_notab(s, c, x, n); }
static void s_bitwise(nn_ptr s, nn_ptr c, nn_srcptr x, slong n, slong num, double w)
{ (void) num; (void) w; fixed_sin_cos_bitwise_rs(s, c, x, n, 0); }
static void s_dio(nn_ptr s, nn_ptr c, nn_srcptr x, slong n, slong num, double w)
{ _fixed_sin_cos_diophantine_tune(s, c, x, n, num, w); }

static double
tmin_exp(exp_fn f, nn_ptr y, nn_srcptr xs, slong n, slong num, double w, int reps)
{
    double best = 1e300, t0;
    int r;
    slong k;
    for (r = 0; r < reps; r++)
        for (k = 0; k < NX; k++)
        {
            t0 = wall();
            f(y, xs + k * n, n, num, w);
            best = FLINT_MIN(best, wall() - t0);
        }
    return best;
}

static double
tmin_sc(sc_fn f, nn_ptr y1, nn_ptr y2, nn_srcptr xs, slong n, slong num,
    double w, int reps)
{
    double best = 1e300, t0;
    int r;
    slong k;
    for (r = 0; r < reps; r++)
        for (k = 0; k < NX; k++)
        {
            t0 = wall();
            f(y1, y2, xs + k * n, n, num, w);
            best = FLINT_MIN(best, wall() - t0);
        }
    return best;
}

static void
run(slong n, double wfac, int reps, flint_rand_t st)
{
    static const slong nums[] = { 13, 20, 32, 48 };
    slong k, i;
    nn_ptr xs = flint_malloc(NX * n * sizeof(ulong));
    nn_ptr y1 = flint_malloc((n + 2) * sizeof(ulong));
    nn_ptr y2 = flint_malloc((n + 2) * sizeof(ulong));
    double te, ts, t0, pre, warm, w = wfac * FLINT_BITS * n;

    for (k = 0; k < NX; k++)
        flint_mpn_urandomb(xs + k * n, st, FLINT_BITS * n);

    /* warm up before the first timing */
    {
        double t0 = wall();
        while (wall() - t0 < 0.05)
        {
            fixed_exp_notab(y1, xs, n);
            fixed_sin_cos_notab(y1, y2, xs, n);
        }
    }

    te = tmin_exp(e_notab, y1, xs, n, 0, 0, reps);
    ts = tmin_sc(s_notab, y1, y2, xs, n, 0, 0, reps);

    flint_printf("n = %wd limbs (%wd bits): exp_notab %.4f s, sin_cos_notab %.4f s,"
        " sin_cos_notab / exp_notab = %.2f\n", n, FLINT_BITS * n, te, ts, ts / te);

    flint_printf("  exp, in multiples of exp_notab:\n    %-30s %9s %9s\n",
        "method", "precomp", "warm");
    _fixed_exp_logs_clear();
    t0 = wall();
    e_bitwise(y1, xs, n, 0, 0);
    pre = wall() - t0;
    warm = tmin_exp(e_bitwise, y1, xs, n, 0, 0, reps);
    flint_printf("    %-30s %9.2f %9.3f   (r = %d)\n", "bitwise", (pre - warm) / te,
        warm / te, fixed_exp_bitwise_rs_default_r(n));
    for (i = 0; i < 4; i++)
    {
        char name[64];
        fixed_rel_table(0, nums[i]);   /* the relation table is not timed */
        _fixed_log_primes_clear();
        t0 = wall();
        e_dio(y1, xs, n, nums[i], w);
        pre = wall() - t0;
        warm = tmin_exp(e_dio, y1, xs, n, nums[i], w, reps);
        flint_sprintf(name, "diophantine, %wd primes, w=%g", nums[i], wfac);
        flint_printf("    %-30s %9.2f %9.3f\n", name, (pre - warm) / te, warm / te);
        fflush(stdout);
    }

    flint_printf("  sin_cos, in multiples of exp_notab (sin_cos_notab = %.2f):\n"
        "    %-30s %9s %9s\n", ts / te, "method", "precomp", "warm");
    _fixed_atans_clear();
    t0 = wall();
    s_bitwise(y1, y2, xs, n, 0, 0);
    pre = wall() - t0;
    warm = tmin_sc(s_bitwise, y1, y2, xs, n, 0, 0, reps);
    flint_printf("    %-30s %9.2f %9.3f   (r = %d)\n", "bitwise", (pre - warm) / te,
        warm / te, fixed_trig_bitwise_rs_default_r(n));
    for (i = 0; i < 4; i++)
    {
        char name[64];
        fixed_rel_table(1, nums[i]);
        _fixed_atan_gauss_clear();
        t0 = wall();
        s_dio(y1, y2, xs, n, nums[i], w);
        pre = wall() - t0;
        warm = tmin_sc(s_dio, y1, y2, xs, n, nums[i], w, reps);
        flint_sprintf(name, "diophantine, %wd primes, w=%g", nums[i], wfac);
        flint_printf("    %-30s %9.2f %9.3f\n", name, (pre - warm) / te, warm / te);
        fflush(stdout);
    }
    flint_printf("\n");

    flint_free(xs);
    flint_free(y1);
    flint_free(y2);
}

int
main(int argc, char ** argv)
{
    flint_rand_t st;
    double wfac = 4.0;
    int reps = 2, i, any = 0;
    static const slong defaults[] = { 256, 1024, 4096, 16384, 65536 };

    flint_rand_init(st);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-w") && i + 1 < argc)
            wfac = atof(argv[++i]);
        else if (!strcmp(argv[i], "-r") && i + 1 < argc)
            reps = atoi(argv[++i]);
    }
    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-w") || !strcmp(argv[i], "-r"))
        {
            i++;
            continue;
        }
        run(atol(argv[i]), wfac, reps, st);
        any = 1;
    }
    if (!any)
        for (i = 0; i < 5; i++)
            run(defaults[i], wfac, reps, st);

    flint_rand_clear(st);
    flint_cleanup_master();
    return 0;
}
