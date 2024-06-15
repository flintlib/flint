/*
    Copyright 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"

#if FLINT_USES_BLAS
#include <cblas.h>
#include "longlong.h"  // for FLINT_BIT_COUNT
#include "fmpz_mat.h"
#include "profiler.h"

int main(void)
{
    slong dim, i, reps;
    flint_bitcnt_t Abits, Bbits, Cbits, maxAbits;
    slong time1, time2, Cbitscheckup = 0;
    double ratio;
    timeit_t timer;
    FLINT_TEST_INIT(state);

    flint_set_num_threads(8);
    //openblas_set_num_threads(8);

    for (dim = 50; dim <= 3000; dim += 2 + dim/4)
    {
        fmpz_mat_t A, B, C, D;

        fmpz_mat_init(A, dim, dim);
        fmpz_mat_init(B, dim, dim);
        fmpz_mat_init(C, dim, dim);
        fmpz_mat_init(D, dim, dim);

        maxAbits = 700 + (20+50)*(3000-800)/(20+dim);

        if (--Cbitscheckup < 0)
        {
            flint_printf("  C bits : ");
            for (Abits = 10; Abits < maxAbits; Abits += 1 + Abits/4)
            {
                Bbits = Abits;
                Cbits = Abits + Bbits + FLINT_BIT_COUNT(dim);
                flint_printf("%5wd  ", Cbits);
            }
            flint_printf("\n");
            Cbitscheckup = 5;
        }

        flint_printf("dim %4wd: |", dim);
        fflush(stdout);

        for (Abits = 10; Abits < maxAbits; Abits += 1 + Abits/4)
        {
            Bbits = Abits;
            Cbits = Abits + Bbits + FLINT_BIT_COUNT(dim);

            fmpz_mat_randtest(A, state, Abits);
            fmpz_mat_randtest(B, state, Bbits);

            reps = 1 + 2000000000/dim/dim/dim/(20 + Cbits);

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul_blas(C, A, B);
            timeit_stop(timer);
            time1 = timer->wall;

            timeit_start(timer);
            for (i = reps; i > 0; i--)
                fmpz_mat_mul(D, A, B);
            timeit_stop(timer);
            time2 = timer->wall;

            ratio = (double)time2/(double)time1;
            ratio = FLINT_MIN(ratio, 9.0);

            flint_printf(" %.2f |", ratio);
            fflush(stdout);
        }

        flint_printf("\n");

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    return 0;
}

#else
int main(void)
{
    return 0;
}
#endif
