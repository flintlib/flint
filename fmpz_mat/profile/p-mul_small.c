/*
    Copyright 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz_mat.h"
#include "profiler.h"

int main(void)
{
    slong m, k, n;
    slong dim, i, reps;
    slong new_total, new_mint, new_maxt;
    slong old_total, old_mint, old_maxt;
    double total_den, den;
    flint_bitcnt_t Abits, Bbits;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(1);

    flint_printf("(m, k, n) *** timings for mul_small nanoseconds per m*k*n (mul in parentheses) ***\n");

    for (m = 3; m < 10; m++)
    for (k = 3; k < 10; k++)
    for (n = 3; n < 10; n++)
    {
        fmpz_mat_t A, B, C, D;

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);
        fmpz_mat_init(D, m, n);

        reps = 1 + 8000000/m/n/k;

        den = reps*m*n*k;

        total_den = 0;

        new_total = 0;
        new_mint = 10000000000;
        new_maxt = 0;

        old_total = 0;
        old_mint = 10000000000;
        old_maxt = 0;

        for (Abits = 14; Abits <= SMALL_FMPZ_BITCOUNT_MAX; Abits += 16)
        for (Bbits = Abits; Bbits <= SMALL_FMPZ_BITCOUNT_MAX; Bbits += 16)
        {
            fmpz_mat_randtest(A, state, Abits);
            fmpz_mat_randtest(B, state, Bbits);

            total_den += den;

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                _fmpz_mat_mul_small(C, A, B);
            timeit_stop(timer);

            new_total += timer->wall;
            new_mint = FLINT_MIN(new_mint, timer->wall);
            new_maxt = FLINT_MAX(new_maxt, timer->wall);

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul(D, A, B);
            timeit_stop(timer);

            old_total += timer->wall;
            old_mint = FLINT_MIN(old_mint, timer->wall);
            old_maxt = FLINT_MAX(old_maxt, timer->wall);
        }

        flint_printf("(%2wd,%2wd,%2wd): min %.3f ns (%.3f)  max %.3f ns (%.3f)  aver %.3f ns (%.3f)  ratio %.3f\n",
                     m, k, n,
                     1000000*new_mint/den, 1000000*old_mint/den,
                     1000000*new_maxt/den, 1000000*old_maxt/den,
                     1000000*new_total/total_den, 1000000*old_total/total_den,
                     (double)new_total/old_total);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    flint_printf("*** timings for mul_small nanoseconds per dim^3 (mul in parentheses) ***\n");

    for (dim = 5; dim <= 1100; dim += 2 + dim/4)
    {
        fmpz_mat_t A, B, C, D;

        fmpz_mat_init(A, dim, dim);
        fmpz_mat_init(B, dim, dim);
        fmpz_mat_init(C, dim, dim);
        fmpz_mat_init(D, dim, dim);

        reps = 1 + 5000000/dim/dim/dim;

        den = reps*dim*dim*dim;

        total_den = 0;

        new_total = 0;
        new_mint = 10000000000;
        new_maxt = 0;

        old_total = 0;
        old_mint = 10000000000;
        old_maxt = 0;

        for (Abits = 14; Abits <= SMALL_FMPZ_BITCOUNT_MAX; Abits += 8)
        for (Bbits = Abits; Bbits <= SMALL_FMPZ_BITCOUNT_MAX; Bbits += 8)
        {
            fmpz_mat_randtest(A, state, Abits);
            fmpz_mat_randtest(B, state, Bbits);

            total_den += den;

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                _fmpz_mat_mul_small(C, A, B);
            timeit_stop(timer);

            new_total += timer->wall;
            new_mint = FLINT_MIN(new_mint, timer->wall);
            new_maxt = FLINT_MAX(new_maxt, timer->wall);

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul(D, A, B);
            timeit_stop(timer);

            old_total += timer->wall;
            old_mint = FLINT_MIN(old_mint, timer->wall);
            old_maxt = FLINT_MAX(old_maxt, timer->wall);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("oops C != D %wu %wu\n", Abits, Bbits);
                flint_abort();
            }
        }

        flint_printf("%4wd: min %.3f ns (%.3f)  max %.3f ns (%.3f)  aver %.3f ns (%.3f)  ratio %.3f\n",
                     dim,
                     1000000*new_mint/den, 1000000*old_mint/den,
                     1000000*new_maxt/den, 1000000*old_maxt/den,
                     1000000*new_total/total_den, 1000000*old_total/total_den,
                     (double)new_total/old_total);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
