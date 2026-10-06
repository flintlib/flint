/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Profile the elementary functions of mp_real against arb across
   thread counts:

       p-threads [FUNCS] [MINBITS] [MAXBITS] [THREADS] [SECONDS]

   FUNCS    comma-separated list among exp, sin_cos, log, atan, or all
            (the default); or "reduced", which instead times the
            algorithms of _mp_real_exp_reduced and
            _mp_real_sin_cos_reduced (0 tuned, 2 series, 3 one burst
            step, 4 full burst) at reduction depths r = 32, 128 and 512
            for each thread count (a sixth argument gives other
            depths, comma-separated), the data for retuning the burst
            thresholds with threads
   MINBITS  smallest precision, default 10^5
   MAXBITS  largest precision, default 10^7 (sizes step by factors of
            about 3.2, i.e. two per decade)
   THREADS  comma-separated thread counts, default 1,2,4,8 (each count
            is used for both libraries through flint_set_num_threads;
            counts beyond the machine's cores only add noise)
   SECONDS  the minimum time spent per measurement, default 1
            (each figure is the best of three such runs)

   For each function, precision and thread count the program prints
   the time per call of arb and of mp_real, the ratio arb / mp_real,
   each library's speedup over its own one-thread time at that
   precision, and (on Linux) the peak memory of each library's calls
   in MB above the level before them -- the resident high-water mark
   of the process, reset before each measurement, which includes the
   FFT's per-thread tables the first call on each thread leaves behind.  Both libraries are given the same argument, a random
   number in [1/4, 1) with a tiny radius, and the reference constants
   (log 2, pi/4 and the tables of the diophantine reductions) are
   built by a first untimed call at the largest precision, so the
   figures are for the evaluation alone.  Results are checked to
   agree between the two libraries. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include "profiler.h"
#include "flint.h"
#include "arb.h"
#include "mpn_extras.h"
#include "mp_real.h"

#define NFUNCS 4

/* peak resident memory in MB (Linux: VmHWM of /proc/self/status, the
   high-water mark reset by writing 5 to /proc/self/clear_refs);
   -1 elsewhere */
static void
peak_reset(void)
{
#ifdef __GLIBC__
    /* return the freed heap to the system, so that the next call's
       footprint shows as fresh pages */
    malloc_trim(0);
#endif
#ifdef __linux__
    FILE * f = fopen("/proc/self/clear_refs", "w");
    if (f != NULL)
    {
        fputs("5\n", f);
        fclose(f);
    }
#endif
}

static double
peak_mb(void)
{
#ifdef __linux__
    FILE * f = fopen("/proc/self/status", "r");
    char line[256];
    long kb = -1;
    if (f == NULL)
        return -1.0;
    while (fgets(line, sizeof(line), f) != NULL)
        if (strncmp(line, "VmHWM:", 6) == 0)
            sscanf(line + 6, "%ld", &kb);
    fclose(f);
    return (kb < 0) ? -1.0 : kb / 1024.0;
#else
    return -1.0;
#endif
}

static double
current_mb(void)
{
#ifdef __linux__
    FILE * f = fopen("/proc/self/status", "r");
    char line[256];
    long kb = -1;
    if (f == NULL)
        return -1.0;
    while (fgets(line, sizeof(line), f) != NULL)
        if (strncmp(line, "VmRSS:", 6) == 0)
            sscanf(line + 6, "%ld", &kb);
    fclose(f);
    return (kb < 0) ? -1.0 : kb / 1024.0;
#else
    return -1.0;
#endif
}

static const char * func_names[NFUNCS] = { "exp", "sin_cos", "log", "atan" };

/* one call of function fi at prec bits (arb or mp_real) */
static void
call_arb(int fi, arb_t r1, arb_t r2, const arb_t x, slong prec)
{
    switch (fi)
    {
        case 0: arb_exp(r1, x, prec); break;
        case 1: arb_sin_cos(r1, r2, x, prec); break;
        case 2: arb_log(r1, x, prec); break;
        default: arb_atan(r1, x, prec); break;
    }
}

static void
call_mp(int fi, mp_real_t r1, mp_real_t r2, const mp_real_t x, slong prec)
{
    switch (fi)
    {
        case 0: mp_real_exp_bits(r1, x, prec); break;
        case 1: mp_real_sin_cos_bits(r1, r2, x, prec); break;
        case 2: mp_real_log_bits(r1, x, prec); break;
        default: mp_real_atan_bits(r1, x, prec); break;
    }
}

/* the best of three runs of at least min_ms milliseconds, in seconds
   per call */
static double
time_calls(int which, int fi, arb_t a1, arb_t a2, const arb_t xa,
    mp_real_t m1, mp_real_t m2, const mp_real_t xm, slong prec, double min_ms)
{
    double best = 1e300;
    int pass;

    for (pass = 0; pass < 3; pass++)
    {
        slong reps = 1;
        timeit_t tm;

        for (;;)
        {
            slong it;
            timeit_start(tm);
            for (it = 0; it < reps; it++)
            {
                if (which == 0)
                    call_arb(fi, a1, a2, xa, prec);
                else
                    call_mp(fi, m1, m2, xm, prec);
            }
            timeit_stop(tm);
            if (tm->wall >= min_ms)
                break;
            reps *= 2;
        }
        best = FLINT_MIN(best, 1e-3 * tm->wall / reps);
    }
    return best;
}

static int
agree(const arb_t a, const mp_real_t m)
{
    arb_t t;
    int ok;
    arb_init(t);
    mp_real_get_arb(t, m);
    ok = arb_overlaps(a, t);
    arb_clear(t);
    return ok;
}

/* the reduced-function algorithms: time per call in seconds */
static double
time_reduced(int trig, nn_ptr y1, nn_ptr y2, nn_srcptr t, slong n,
    flint_bitcnt_t r, int alg, double min_ms)
{
    double best = 1e300;
    int pass;

    for (pass = 0; pass < 3; pass++)
    {
        slong reps = 1;
        timeit_t tm;

        for (;;)
        {
            slong it;
            timeit_start(tm);
            for (it = 0; it < reps; it++)
            {
                if (trig)
                    _mp_real_sin_cos_reduced(y1, y2, NULL, t, n, r, alg);
                else
                    _mp_real_exp_reduced(y1, NULL, t, n, r, alg);
            }
            timeit_stop(tm);
            if (tm->wall >= min_ms)
                break;
            reps *= 2;
        }
        best = FLINT_MIN(best, 1e-3 * tm->wall / reps);
    }
    return best;
}

static void
profile_reduced(slong minbits, slong maxbits, const slong * threads,
    slong nthreads, double min_ms, const slong * rs, slong nrs,
    flint_rand_t state)
{
    static const int algs[] = { 0, 2, 3, 4 };
    slong bits, ri, ti, i;
    int trig, ai;

    flint_printf("%8s %8s %5s %7s %12s %12s %12s %12s\n", "func", "bits", "r",
        "threads", "alg 0 (s)", "alg 2 (s)", "alg 3 (s)", "alg 4 (s)");

    for (trig = 0; trig < 2; trig++)
    for (bits = minbits; bits <= maxbits; bits = (slong) (bits * 3.162278 + 1))
    for (ri = 0; ri < nrs; ri++)
    {
        slong n = bits / FLINT_BITS + 1;
        flint_bitcnt_t r = (flint_bitcnt_t) rs[ri];
        nn_ptr t, y1, y2;

        t = flint_malloc(n * sizeof(ulong));
        y1 = flint_malloc((n + 2) * sizeof(ulong));
        y2 = flint_malloc((n + 2) * sizeof(ulong));
        for (i = 0; i < n; i++)
            t[i] = n_randtest(state);
        /* t < 2^-r */
        for (i = 0; i < (slong) (r / FLINT_BITS) && i < n; i++)
            t[n - 1 - i] = 0;
        if ((r % FLINT_BITS) && (slong) (r / FLINT_BITS) < n)
            t[n - 1 - r / FLINT_BITS] >>= (r % FLINT_BITS);

        for (ti = 0; ti < nthreads; ti++)
        {
            flint_set_num_threads((int) threads[ti]);
            flint_printf("%8s %8wd %5wd %7wd", trig ? "sin_cos" : "exp",
                bits, (slong) r, threads[ti]);
            for (ai = 0; ai < 4; ai++)
                flint_printf(" %12.5f", time_reduced(trig, y1, y2, t, n, r,
                    algs[ai], min_ms));
            flint_printf("\n");
            fflush(stdout);
        }

        flint_free(t);
        flint_free(y1);
        flint_free(y2);
    }
}

int
main(int argc, char * argv[])
{
    int use[NFUNCS] = { 1, 1, 1, 1 }, reduced = 0;
    slong minbits = 100000, maxbits = 10000000;
    slong threads[64], nthreads = 0, rs[64], nrs = 0;
    double min_ms = 1000.0;
    slong bits, i, fi, ti;
    flint_rand_t state;

    if (argc > 1 && strcmp(argv[1], "reduced") == 0)
        reduced = 1;
    else if (argc > 1 && strcmp(argv[1], "all") != 0)
    {
        for (fi = 0; fi < NFUNCS; fi++)
            use[fi] = 0;
        for (fi = 0; fi < NFUNCS; fi++)
        {
            const char * p = strstr(argv[1], func_names[fi]);
            /* "sin_cos" contains no other name; "exp" etc. are distinct */
            if (p != NULL)
                use[fi] = 1;
        }
    }
    if (argc > 2) minbits = atol(argv[2]);
    if (argc > 3) maxbits = atol(argv[3]);
    if (argc > 4)
    {
        char * s = argv[4];
        while (*s && nthreads < 64)
        {
            threads[nthreads++] = atol(s);
            while (*s && *s != ',') s++;
            if (*s == ',') s++;
        }
    }
    else
    {
        threads[0] = 1; threads[1] = 2; threads[2] = 4; threads[3] = 8;
        nthreads = 4;
    }
    if (argc > 5) min_ms = 1000.0 * atof(argv[5]);
    if (argc > 6)
    {
        char * s = argv[6];
        while (*s && nrs < 64)
        {
            rs[nrs++] = atol(s);
            while (*s && *s != ',') s++;
            if (*s == ',') s++;
        }
    }
    else
    {
        rs[0] = 32; rs[1] = 128; rs[2] = 512;
        nrs = 3;
    }

    flint_rand_init(state);

    if (reduced)
    {
        profile_reduced(minbits, maxbits, threads, nthreads, min_ms, rs, nrs,
            state);
        flint_set_num_threads(1);
        flint_rand_clear(state);
        flint_cleanup_master();
        return 0;
    }

    flint_printf("%8s %10s %7s %12s %12s %8s %7s %7s %8s %8s\n", "func", "bits",
        "threads", "arb (s)", "mp_real (s)", "arb/mp", "arb x", "mp x",
        "arb MB", "mp MB");

    for (fi = 0; fi < NFUNCS; fi++)
    {
        if (!use[fi])
            continue;

        for (bits = minbits; bits <= maxbits; bits = (slong) (bits * 3.162278 + 1))
        {
            slong n = bits / FLINT_BITS + 1;
            arb_t xa, a1, a2;
            mp_real_t xm, m1, m2;
            nn_ptr v;
            fmpz_t f;
            double ta1 = 0.0, tm1 = 0.0;

            /* x in [1/4, 1) with a radius of 2^-(bits+10), the same
               ball for both libraries */
            v = flint_malloc(n * sizeof(ulong));
            for (i = 0; i < n; i++)
                v[i] = n_randtest(state);
            v[n - 1] |= (UWORD(1) << (FLINT_BITS - 2));
            fmpz_init(f);
            fmpz_set_ui_array(f, v, n);
            arb_init(xa); arb_init(a1); arb_init(a2);
            arb_set_fmpz(xa, f);
            arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
            arb_add_error_2exp_si(xa, -bits - 10);
            mp_real_init(xm); mp_real_init(m1); mp_real_init(m2);
            _mp_real_set_mpn_2exp(xm, v, n, -FLINT_BITS * n);
            mp_real_add_error_2exp_si(xm, -bits - 10);

            for (ti = 0; ti < nthreads; ti++)
            {
                double ta, tm, ma, mm, base;

                flint_set_num_threads((int) threads[ti]);

                /* warm up: constants and tables, outside the timing */
                call_arb(fi, a1, a2, xa, bits);
                call_mp(fi, m1, m2, xm, bits);
                if (!agree(a1, m1) || (fi == 1 && !agree(a2, m2)))
                    flint_printf("WARNING: %s at %wd bits, %wd threads: "
                        "results disagree\n", func_names[fi], bits, threads[ti]);

                peak_reset();
                base = current_mb();
                ta = time_calls(0, fi, a1, a2, xa, m1, m2, xm, bits, min_ms);
                ma = peak_mb() - base;
                peak_reset();
                base = current_mb();
                tm = time_calls(1, fi, a1, a2, xa, m1, m2, xm, bits, min_ms);
                mm = peak_mb() - base;
                if (ti == 0)
                {
                    ta1 = ta;
                    tm1 = tm;
                }

                flint_printf("%8s %10wd %7wd %12.4f %12.4f %8.2f %7.2f %7.2f %8.0f %8.0f\n",
                    func_names[fi], bits, threads[ti], ta, tm, ta / tm,
                    ta1 / ta, tm1 / tm, ma, mm);
                fflush(stdout);
            }

            arb_clear(xa); arb_clear(a1); arb_clear(a2);
            mp_real_clear(xm); mp_real_clear(m1); mp_real_clear(m2);
            fmpz_clear(f);
            flint_free(v);
        }
    }

    flint_set_num_threads(1);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
