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
#include "profiler.h"
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "test_helpers.h"


int main(void)
{
    slong i, j, reps;
    slong dim, sz;
    fmpz_mat_t A, B, C, D;
    fmpz_t t;
    timeit_t timer;
    slong time1, time2;

    reps = 1;

    for (dim = 6; dim > 0; dim--)
    {
        flint_printf("****** %wd x %wd *********\n", dim, dim);
        for (sz = 10000; sz < 100000; sz += 10000)
        {
            fmpz_mat_init(A, dim, dim);
            fmpz_mat_init(B, dim, dim);
            fmpz_mat_init(C, dim, dim);
            fmpz_mat_init(D, dim, dim);
            fmpz_init_set_ui(t, 1);

            flint_printf("sz = %wd\n", sz);

            fmpz_mul_2exp(t, t, FLINT_BITS*sz);
            fmpz_sub_ui(t, t, 1);

            for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
            {
                fmpz_set(fmpz_mat_entry(A, i, j), t);
                fmpz_set(fmpz_mat_entry(B, i, j), t);
            }

            for (i = 0; i < reps; i++)
            {
                timeit_start(timer);
                fmpz_mat_mul_fft(C, A, B);
                timeit_stop(timer);
                time1 = timer->wall;

                timeit_start(timer);
                fmpz_mat_mul_classical(D, A, B);
                timeit_stop(timer);
                time2 = timer->wall;

                flint_printf("new: %wd,  old: %wd,  ratio: %f\n", time1, time2, (double)time2/(double)time1);

                FLINT_TEST(fmpz_mat_equal(C, D));
            }

            fmpz_clear(t);
            fmpz_mat_clear(A);
            fmpz_mat_clear(B);
            fmpz_mat_clear(C);
            fmpz_mat_clear(D);
        }
    }

    return 0;
}

