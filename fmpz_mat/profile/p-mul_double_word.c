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
    slong sign;
    slong dim, i, reps, total, mint, maxt;
    double total_den, den;
    flint_bitcnt_t Abits, Bbits;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(1);

    flint_printf("*** timings are nanoseconds per dim^3 ***\n");

    for (sign = 0; sign <= 1; sign++)
    for (dim = 5; dim <= 400; dim += 4 + dim/16)
    {
        fmpz_mat_t A, B, C, D, E;

        fmpz_mat_init(A, dim, dim);
        fmpz_mat_init(B, dim, dim);
        fmpz_mat_init(C, dim, dim);
        fmpz_mat_init(D, dim, dim);
        fmpz_mat_init(E, dim, dim);

        reps = 2 + 3000000/dim/dim/dim;

        den = reps*dim*dim*dim;

        total = total_den = 0;
        mint = 10000000000;
        maxt = 0;

        for (Abits = FLINT_BITS - 1; Abits + sign <= 2*FLINT_BITS; Abits += 1)
        for (Bbits = Abits; Bbits + sign <= 2*FLINT_BITS; Bbits += 27)
        {
            if (sign)
            {
                fmpz_mat_randtest(A, state, Abits);
                fmpz_mat_randtest(B, state, Bbits);
            }
            else
            {
                fmpz_mat_randtest_unsigned(A, state, Abits);
                fmpz_mat_randtest_unsigned(B, state, Bbits);                
            }

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                _fmpz_mat_mul_double_word(E, A, B);
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
                    flint_printf("E: ");
                    fmpz_mat_print_pretty(E);
                    flint_printf("\n");

                    flint_printf("D: ");
                    fmpz_mat_print_pretty(D);
                    flint_printf("\n");

                    flint_printf("oops %d %wu %wu\n", sign, Abits, Bbits);
                    flint_abort();
                }
            }
        }

        flint_printf("sign %d dim %3wd: min %.3f ns  max %.3f ns  aver %.3f ns\n",
             sign, dim, 1000000*mint/den, 1000000*maxt/den, 1000000*total/total_den);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(E);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}
