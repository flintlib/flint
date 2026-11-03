/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Check the tuning of Khinchin's constant (_khinchin_params in
   const_khinchin.c) against its neighbours:

       p-const_khinchin_tune [MINBITS] [MAXBITS] [THREADS] [SECONDS]

   MINBITS  smallest precision, default 1000
   MAXBITS  largest precision, default 300000 (sizes step by factors of
            about 3.2, i.e. two per decade)
   THREADS  comma-separated thread counts, default 1,2,4,8 (counts
            beyond the machine's cores only add noise)
   SECONDS  the minimum time spent per measurement, default 0.3 (each
            figure is the best of three such runs)

   The parameters are x = log2 K(n1), K(n1) the largest k summed
   directly (the Bernoulli columns being n < n1 = F / (2 x) for F
   fraction bits), K0, the logarithm rows being k < K0, and nb, the
   number of Bernoulli chunks run at once (each with about THREADS / nb
   threads inside). For each precision and thread count the program
   prints the time with the default parameters, then the time relative
   to it with x scaled by 0.7 ... 1.4, with K0 from 4 to 24, and with
   nb = 1, 2, 4, ... up to the thread count (the other parameters at
   their defaults): entries below 1 point to a better choice. The
   default x is recovered from the default n1, so it is rounded. The
   memory of the Bernoulli chunks grows with nb; the default nb is the
   largest whose power tables stay below about 1 GiB. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "profiler.h"
#include "flint.h"
#include "thread_support.h"
#include "mp_real.h"

/* the best of three runs of at least min_ms milliseconds, in seconds
   per call */
static double
time_calls(slong n, slong K0, slong n1, slong nb, double min_ms)
{
    double best = 1e300;
    int pass;
    mp_real_t res;

    mp_real_init(res);
    for (pass = 0; pass < 3; pass++)
    {
        slong reps = 1, it;
        timeit_t tm;

        for (;;)
        {
            timeit_start(tm);
            for (it = 0; it < reps; it++)
                _mp_real_const_khinchin_tune(res, n, K0, n1, nb);
            timeit_stop(tm);
            if (tm->wall >= min_ms)
                break;
            reps *= 2;
        }
        best = FLINT_MIN(best, 1e-3 * tm->wall / reps);
    }
    mp_real_clear(res);
    return best;
}

static slong
parse_list(slong * out, slong maxlen, const char * s)
{
    slong len = 0;
    while (*s != '\0' && len < maxlen)
    {
        out[len++] = strtol(s, (char **) &s, 10);
        if (*s == ',')
            s++;
    }
    return len;
}

int
main(int argc, char ** argv)
{
    slong minbits = 1000, maxbits = 300000, threads[16] = {1, 2, 4, 8};
    slong nthreads = 4, ti, bits, i;
    double min_ms = 300.0;
    static const double xs[] = {0.7, 0.8, 0.9, 1.1, 1.2, 1.4};
    static const slong ks[] = {4, 6, 8, 10, 12, 16, 24};

    if (argc > 1) minbits = atol(argv[1]);
    if (argc > 2) maxbits = atol(argv[2]);
    if (argc > 3) nthreads = parse_list(threads, 16, argv[3]);
    if (argc > 4) min_ms = 1000.0 * atof(argv[4]);

    for (bits = minbits; bits <= maxbits; bits = (slong) (bits * 3.1623 + 0.5))
    {
        slong n = (bits + FLINT_BITS - 1) / FLINT_BITS;
        double F = (double) FLINT_BITS * n;

        for (ti = 0; ti < nthreads; ti++)
        {
            slong K0, n1, nb;
            double x, t0;

            flint_set_num_threads(threads[ti]);
            _mp_real_const_khinchin_params(&K0, &n1, &nb, n);
            x = F / (2.0 * n1);
            time_calls(n, K0, n1, nb, 0.0);     /* warm up */
            t0 = time_calls(n, K0, n1, nb, min_ms);

            flint_printf("bits %wd, threads %wd: K0 = %wd, x = %.2f (n1 = %wd), nb = %wd: %.4g s\n",
                bits, threads[ti], K0, x, n1, nb, t0);

            flint_printf("    x * ");
            for (i = 0; i < (slong) (sizeof(xs) / sizeof(xs[0])); i++)
            {
                slong m = FLINT_MAX((slong) (F / (2.0 * x * xs[i])), 1);
                flint_printf("  %.1f: %.3f", xs[i], time_calls(n, K0, m, nb, min_ms) / t0);
            }
            flint_printf("\n    K0 =");
            for (i = 0; i < (slong) (sizeof(ks) / sizeof(ks[0])); i++)
                flint_printf("  %wd: %.3f", ks[i], time_calls(n, ks[i], n1, nb, min_ms) / t0);
            if (threads[ti] > 1)
            {
                slong b;
                flint_printf("\n    nb =");
                for (b = 1; b < threads[ti]; b *= 2)
                    flint_printf("  %wd: %.3f", b, time_calls(n, K0, n1, b, min_ms) / t0);
                flint_printf("  %wd: %.3f", threads[ti], time_calls(n, K0, n1, threads[ti], min_ms) / t0);
            }
            flint_printf("\n");
        }
    }

    flint_set_num_threads(1);
    flint_cleanup_master();
    return 0;
}
