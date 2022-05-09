/*
    Copyright 2009 William Hart
    Copyright 2010,2011 Fredrik Johansson
    Copyright 2014 Abhinav Baid

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
#include "fmpz_lll.h"
#include "fmpz.h"
#include "ulong_extras.h"

typedef struct
{
    slong dim;
    int algorithm;
} mat_lll_t;


void
sample(void *arg, ulong count)
{
    mat_lll_t *params = (mat_lll_t *) arg;
    slong i, dim = params->dim;
    int algorithm = params->algorithm;
    fmpq_t delta, eta;
    fmpz_lll_t fl;

    fmpz_mat_t A, B, C, D;
    FLINT_TEST_INIT(state);


    fmpz_mat_init(A, dim, dim);
    fmpq_init(delta);
    fmpq_init(eta);

    fmpq_set_si(delta, 3, 4);
    fmpq_set_si(eta, 81, 100);
    fmpz_lll_context_init(fl, 0.75, 0.81, 1, 0);

    fmpz_mat_randajtai(A, state, 0.5);
    fmpz_mat_init_set(B, A);
    fmpz_mat_init_set(C, A);
    fmpz_mat_init_set(D, A);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
        {
            fmpz_mat_lll_original(A, delta, eta);
        }
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
        {
            fmpz_mat_lll_storjohann(B, delta, eta);
        }
    else if (algorithm == 2)
        for (i = 0; i < count; i++)
        {
            fmpz_lll_wrapper(C, NULL, fl);
        }
    else if (algorithm == 3)
        for (i = 0; i < count; i++)
        {
            fmpz_lll(D, NULL, fl);
        }

    prof_stop();

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
    fmpz_mat_clear(D);
    fmpq_clear(delta);
    fmpq_clear(eta);
    flint_randclear(state);
}

int
main(void)
{
    double min_classical, min_storjohann, min_wrapper, min_default, max;
    mat_lll_t params;
    slong dim;

    flint_printf("fmpz_lll :\n");

    for (dim = 100; dim <= 500; dim += 10)
    {
        params.dim = dim;

        params.algorithm = 0;
        prof_repeat(&min_classical, &max, sample, &params);

        params.algorithm = 1;
        prof_repeat(&min_storjohann, &max, sample, &params);

        params.algorithm = 2;
        prof_repeat(&min_wrapper, &max, sample, &params);

        params.algorithm = 3;
        prof_repeat(&min_default, &max, sample, &params);

        flint_printf
            ("dim = %wd classical/storjohann/wrapper/default %.2f %.2f %.2f %.2f (us)\n",
             dim, min_classical, min_storjohann, min_wrapper, min_default);
    }

    return 0;
}
