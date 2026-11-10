/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Compare the two algorithms for Khinchin's constant, for choosing
   MP_REAL_KHINCHIN_AFE_CUTOFF (the zeta-value sum below it, the theta
   splitting from it):

       p-const_khinchin_afe [DIGITS] [THREADS] [A] [ZETA]

   DIGITS   comma-separated decimal precisions, default
            10000,30000,100000,300000
   THREADS  comma-separated thread counts, default 1
   A        comma-separated split points u0 = 4^-a of the theta
            splitting, default 2 (the default of _mp_real_const_khinchin_afe)
   ZETA     0 to skip the zeta-value sum (for sizes where it is known to
            be slower), default 1

   Each figure is the wall time of one run, in seconds; the last column
   is the time of the zeta-value sum divided by that of the theta
   splitting with the first split point (above 1: the theta splitting is
   faster). The results are checked for overlap. Above 10^6 digits
   the runs take hours and the memory of the zeta-value sum becomes
   large; set FLINT_MP_REAL_VERBOSE=1 for progress reports. */

#include <stdlib.h>
#include <stdio.h>
#include "profiler.h"
#include "flint.h"
#include "thread_support.h"
#include "arb.h"
#include "mp_real.h"

static int
_overlaps(const mp_real_t x, const mp_real_t y)
{
    arb_t a, b;
    int r;
    arb_init(a);
    arb_init(b);
    mp_real_get_arb(a, x);
    mp_real_get_arb(b, y);
    r = arb_overlaps(a, b);
    arb_clear(a);
    arb_clear(b);
    return r;
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
    slong digits[32] = {10000, 30000, 100000, 300000};
    slong threads[16] = {1}, as[8] = {2};
    slong nd = 4, nth = 1, na = 1, di, ti, ai;
    int zeta = 1;

    if (argc > 1) nd = parse_list(digits, 32, argv[1]);
    if (argc > 2) nth = parse_list(threads, 16, argv[2]);
    if (argc > 3) na = parse_list(as, 8, argv[3]);
    if (argc > 4) zeta = atoi(argv[4]);

    for (di = 0; di < nd; di++)
    {
        /* digits log2(10) bits */
        slong n = (slong) (digits[di] * 3.3219280948873622 / FLINT_BITS) + 2;

        for (ti = 0; ti < nth; ti++)
        {
            mp_real_t x, y;
            timeit_t tm;
            double tz = 0.0, t1 = 0.0;

            flint_set_num_threads(threads[ti]);
            mp_real_init(x);
            mp_real_init(y);

            flint_printf("digits %wd (%wd limbs), threads %wd:", digits[di], n, threads[ti]);
            fflush(stdout);

            if (zeta)
            {
                timeit_start(tm);
                _mp_real_const_khinchin_tune(y, n, 0, 0, 0);
                timeit_stop(tm);
                tz = 1e-3 * tm->wall;
                flint_printf("  zeta %.3f", tz);
                fflush(stdout);
            }

            for (ai = 0; ai < na; ai++)
            {
                double t;

                timeit_start(tm);
                _mp_real_const_khinchin_afe(x, n, as[ai], 0);
                timeit_stop(tm);
                t = 1e-3 * tm->wall;
                if (ai == 0)
                    t1 = t;
                flint_printf("  afe(a = %wd) %.3f", as[ai], t);
                fflush(stdout);

                /* against the zeta-value sum, else the previous split
                   point */
                if ((zeta || ai > 0) && !_overlaps(x, y))
                {
                    flint_printf("\nFAIL: the results do not overlap\n");
                    flint_abort();
                }
                if (!zeta)
                    mp_real_swap(x, y);
            }

            if (zeta && t1 > 0.0)
                flint_printf("  ratio %.2f", tz / t1);
            flint_printf("\n");

            mp_real_clear(x);
            mp_real_clear(y);
        }
    }

    flint_cleanup_master();
    return 0;
}
