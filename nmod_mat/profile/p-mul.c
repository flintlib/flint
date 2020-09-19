/*
    Copyright 2009 William Hart
    Copyright 2010 Fredrik Johansson

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
#include "nmod_mat.h"
#include "ulong_extras.h"

typedef struct
{
    ulong dim;
    mp_limb_t modulus;
    int algorithm;
} mat_mul_t;

void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    mp_limb_t n = params->modulus;
    ulong i, dim = params->dim;
    int algorithm = params->algorithm;

    nmod_mat_t A, B, C;

    nmod_mat_init(A, dim, dim, n);
    nmod_mat_init(B, dim, dim, n);
    nmod_mat_init(C, dim, dim, n);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
            nmod_mat_mul(C, A, B);
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
            nmod_mat_mul_classical(C, A, B);
    else
        for (i = 0; i < count; i++)
            nmod_mat_mul_strassen(C, A, B);

    prof_stop();

    if (C->entries[0] == WORD(5893479483)) abort();

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

int main(void)
{
    double min_classical, min_strassen, max;
    mat_mul_t params;
    slong dim;
    slong i, k, m, n;
    nmod_mat_t A, B, C;
    flint_rand_t state;
    timeit_t timer1, timer2;

    flint_printf("nmod_mat_mul:\n");

    params.modulus = 40000;

    for (dim = 2; dim <= 512; dim = (slong) ((double) dim * 1.1) + 1)
    {
        params.dim = dim;

        params.algorithm = 1;
        prof_repeat(&min_classical, &max, sample, &params);

        params.algorithm = 2;
        prof_repeat(&min_strassen, &max, sample, &params);

        flint_printf("dim = %wd, classical %.2f us strassen %.2f us\n", 
            dim, min_classical, min_strassen);
    }

    flint_randinit(state);

    for (k = 5; k < 600; k += 2 + k/4)
    {
        flint_printf("k = %wd:\n", k);

        for (m = 5; m < 600; m += 2 + m/4)
        {
            flint_printf("m =%04wd: ", m);

            for (n = 5; n < 600; n += 2 + n/4)
            {
                slong reps = 5 + 4000000/k/m/n;
                nmod_mat_init(A, m, k, 40000);
                nmod_mat_init(B, k, n, 40000);
                nmod_mat_init(C, m, n, 40000);
                nmod_mat_randfull(A, state);
                nmod_mat_randfull(B, state);
                nmod_mat_randfull(C, state);
                timeit_start(timer1);
                for (i = reps; i >= 0; i--)
                    nmod_mat_mul(C, A, B);
                timeit_stop(timer1);
                timeit_start(timer2);
                for (i = reps; i >= 0; i--)
                    nmod_mat_mul_classical(C, A, B);
                timeit_stop(timer2);
                flint_printf(" %0.3f", (double)timer1->wall/(double)FLINT_MAX(1, timer2->wall));
                nmod_mat_clear(A);
                nmod_mat_clear(B);
                nmod_mat_clear(C);
            }

            flint_printf("\n");
        }

        flint_printf("\n");
    }

    flint_randclear(state);

    return 0;
}
