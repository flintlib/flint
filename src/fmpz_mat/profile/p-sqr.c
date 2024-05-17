/*
    Copyright (C) 2007,2015 Anubhav Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "fmpz_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"

typedef struct
{
    slong m;
    slong n;
    slong k;
    int algorithm;
    slong bits;
} mat_mul_t;


void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    slong i, n = params->n, m = params->m ;
    slong bits = params->bits;
    int algorithm = params->algorithm;
    fmpz_mat_t A, B, C;
    FLINT_TEST_INIT(state);


    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, m, n);
    fmpz_mat_init(C, m, n);

    fmpz_mat_randbits(A, state, bits);
    fmpz_mat_set(B, A);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
            fmpz_mat_mul(C, A, B);
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
            fmpz_mat_sqr(C, A);


    prof_stop();

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);

    flint_rand_clear(state);
}

int main(void)
{
    double min_default, min_classical, max;
    mat_mul_t params;
    slong bits, dim;

    for (bits = 32; bits <= 65536; bits *= 2)
    {
        params.bits = bits;

        flint_printf("bits = %wd :\n", params.bits);

        for (dim = 1; dim <= 200; dim = FLINT_MAX(dim * 1.5, dim + 1))
        {
            params.n = dim;
            params.m = dim;

            params.algorithm = 0;
            prof_repeat(&min_default, &max, sample, &params);

            params.algorithm = 1;
            prof_repeat(&min_classical, &max, sample, &params);


            flint_printf("dim = %wd fmpz_mat_mul : %.2f fmpz_mat_sqr : %.2f         ", dim, min_default, min_classical);
            flint_printf("Ratio : %0.2f \n",min_default/min_classical);

        }
    }

    return 0;
}
