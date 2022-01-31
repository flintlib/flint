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
    slong dim, i, reps, den;
    slong min1, max1, total1, min2, max2, total2, count;
    flint_bitcnt_t Abits, Bbits;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(1);

    flint_printf("*** timings are nanoseconds per dim^3 ***\n");

    for (dim = 50; dim <= 5000; dim += 2 + dim/4)
    {
        fmpz_mat_t A, B, C, D;

        fmpz_mat_init(A, dim, dim);
        fmpz_mat_init(B, dim, dim);
        fmpz_mat_init(C, dim, dim);
        fmpz_mat_init(D, dim, dim);

        reps = 1 + 3000000/dim/dim/dim;
        den = reps*dim*dim*dim;
        count = 0;

        total1 = 0;
        min1 = 10000000000;
        max1 = 0;

        total2 = 0;
        min2 = 10000000000;
        max2 = 0;

        for (Abits = FLINT_BITS/2; Abits <= SMALL_FMPZ_BITCOUNT_MAX; Abits += FLINT_BITS)
        for (Bbits = Abits; Bbits <= SMALL_FMPZ_BITCOUNT_MAX; Bbits += 16)
        {
            if (FLINT_BIT_COUNT(dim) + Abits + Bbits <= SMALL_FMPZ_BITCOUNT_MAX)
                continue;

            fmpz_mat_randtest(A, state, Abits);
            fmpz_mat_randtest(B, state, Bbits);

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                _fmpz_mat_mul_small(C, A, B);
            timeit_stop(timer);
            total1 += timer->wall;
            min1 = FLINT_MIN(min1, timer->wall);
            max1 = FLINT_MAX(max1, timer->wall);

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul_multi_mod(D, A, B);
            timeit_stop(timer);
            total2 += timer->wall;
            min2 = FLINT_MIN(min2, timer->wall);
            max2 = FLINT_MAX(max2, timer->wall);

            count++;

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("oops\n");
                flint_abort();
            }
        }

        flint_printf("dim %3wd: "
                     "| min %.3f v %0.3f "
                     "| avr %.3f v %0.3f "
                     "| max %.3f v %0.3f |\n",
                    dim,
                    1000000.0*min1/den, 1000000.0*min2/den,
                    1000000.0*total1/den/count, 1000000.0*total2/den/count,
                    1000000.0*max1/den, 1000000.0*max2/den);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
