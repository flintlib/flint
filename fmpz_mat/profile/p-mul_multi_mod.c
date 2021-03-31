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
    slong dim, i, reps, total, mint, maxt;
    double total_den, den;
    flint_bitcnt_t Abits, Bbits, Cbits;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(4);

    for (dim = 16; dim <= 2000; dim += 5 + dim/8)
    {
        fmpz_mat_t A, B, C, D, E;

        fmpz_mat_init(A, dim, dim);
        fmpz_mat_init(B, dim, dim);
        fmpz_mat_init(C, dim, dim);
        fmpz_mat_init(D, dim, dim);
        fmpz_mat_init(E, dim, dim);

        reps = 1 + 2000000/dim/dim/dim;

        den = reps*dim*dim*dim;

        total = total_den = 0;
        mint = 10000000000;
        maxt = 0;

        for (Abits = 80; Abits < 3*FLINT_BITS; Abits += 45)
        for (Bbits = Abits; Bbits < 3*FLINT_BITS; Bbits += 45)
        {
            Cbits = Abits + Bbits + FLINT_BIT_COUNT(dim) + 1;

            if (Cbits >= 400*FLINT_BITS)
                continue;

            fmpz_mat_randtest(A, state, Abits);
            fmpz_mat_randtest(B, state, Bbits);

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul_multi_mod(E, A, B);
            timeit_stop(timer);

            total += timer->wall;
            total_den += ((double)reps)*dim*dim*dim;
            mint = FLINT_MIN(mint, timer->wall);
            maxt = FLINT_MAX(maxt, timer->wall);

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

        flint_printf("dim %3wd: min %.3f ns  max %.3f ns  aver %.3f ns\n", dim,
                     1000000*mint/den, 1000000*maxt/den, 1000000*total/total_den);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(E);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
