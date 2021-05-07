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
    slong dim, i, reps, total, mint, maxt, count, den;
    flint_bitcnt_t Abits, Bbits;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(8);

    flint_printf("*** timings are nanoseconds per dim^3 ***\n");

    for (dim = 50; dim <= 1000; dim += 5 + dim/8)
    {
        fmpz_mat_t A, B, C, D, E;

        fmpz_mat_init(A, dim, dim);
        fmpz_mat_init(B, dim, dim);
        fmpz_mat_init(C, dim, dim);
        fmpz_mat_init(D, dim, dim);
        fmpz_mat_init(E, dim, dim);

        reps = 1 + 2000000/dim/dim/dim;

        den = reps*dim*dim*dim;
        count = 0;

        total = 0;
        mint = 10000000000;
        maxt = 0;

        for (Abits = 1*FLINT_BITS; Abits < 7*FLINT_BITS; Abits += 40)
        for (Bbits = Abits; Bbits < 7*FLINT_BITS; Bbits += 40)
        {
            fmpz_mat_randtest(A, state, Abits);
            fmpz_mat_randtest(B, state, Bbits);

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul_multi_mod(E, A, B);
            timeit_stop(timer);

            total += timer->wall;
            mint = FLINT_MIN(mint, timer->wall);
            maxt = FLINT_MAX(maxt, timer->wall);
            count++;

            if (dim < 150)
            {
                fmpz_mat_mul_classical_inline(D, A, B);

                if (!fmpz_mat_equal(D, E))
                {
                    flint_printf("oops %wu %wu\n", Abits, Bbits);
                    flint_abort();
                }
            }
        }

        flint_printf("dim %3wd: | min %.3f | aver %.3f | max %.3f ns | %wd reps\n",
                     dim,
                     1000000.0*mint/den,
                     1000000.0*total/den/count,
                     1000000.0*maxt/den,
                     reps);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(E);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
