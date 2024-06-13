/*
    Copyright 2009 William Hart
    Copyright 2010 Fredrik Johansson
    Copyright 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "thread_support.h"

#if FLINT_USES_BLAS
# include <cblas.h>
#endif

typedef struct
{
    slong dim_m;
    slong dim_n;
    slong dim_k;
    ulong modulus;
    int algorithm;
} mat_mul_t;

void sample(void * arg, ulong count)
{
    mat_mul_t * params = (mat_mul_t *) arg;
    int algorithm = params->algorithm;
    nmod_mat_t A, B, C;
    ulong i;
    flint_rand_t state;

    flint_rand_init(state);

    nmod_mat_init(A, params->dim_m, params->dim_k, params->modulus);
    nmod_mat_init(B, params->dim_k, params->dim_n, params->modulus);
    nmod_mat_init(C, params->dim_m, params->dim_n, params->modulus);
    nmod_mat_randfull(A, state);
    nmod_mat_randfull(B, state);
    nmod_mat_randfull(C, state);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
            nmod_mat_mul(C, A, B);
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
            nmod_mat_mul_classical(C, A, B);
    else if (algorithm == 2)
        for (i = 0; i < count; i++)
            nmod_mat_mul_classical_threaded(C, A, B);
    else if (algorithm == 3)
        for (i = 0; i < count; i++)
            nmod_mat_mul_blas(C, A, B);
    else
        for (i = 0; i < count; i++)
            nmod_mat_mul_strassen(C, A, B);

    prof_stop();

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);

    flint_rand_clear(state);
}

int main(void)
{
    double max;
    mat_mul_t params;
    slong dim, i, flint_num, blas_num;

    flint_printf("nmod_mat_mul:\n");

    for (dim = 2; dim <= 100; dim += dim/4 + 1)
    {
        double min_classical, min_strassen;

        params.dim_m = dim;
        params.dim_n = dim;
        params.dim_k = dim;
        params.modulus = 40000;

        params.algorithm = 1;
        prof_repeat(&min_classical, &max, sample, &params);

        params.algorithm = 4;
        prof_repeat(&min_strassen, &max, sample, &params);

        flint_printf("dim = %wd, classical %.2f us strassen %.2f us\n",
                                             dim, min_classical, min_strassen);
    }

    /* output floating point ratios time(mul_blas)/time(mul_blas) */
    for (dim = 200; dim <= 1200; dim += 200)
    {
        flint_printf("dimension %wd\n", dim);

        for (flint_num = 2; flint_num <= 8; flint_num += 1)
        {
            flint_set_num_threads(flint_num);

            for (blas_num = flint_num; blas_num <= flint_num; blas_num *= 2)
            {
                double min_old, min_new, min_ratio = 100;

#if FLINT_USES_BLAS
                //openblas_set_num_threads(blas_num);
#endif

                flint_printf("[flint %wd, blas %wd]: (", flint_num, blas_num);

                for (i = 7; i < FLINT_BITS; i += 8)
                {
                    params.dim_m = dim;
                    params.dim_n = dim;
                    params.dim_k = dim;
                    params.modulus = 2*(UWORD(1) << i) - 1;

                    params.algorithm = 2;
                    prof_repeat(&min_old, &max, sample, &params);

                    params.algorithm = 0;
                    prof_repeat(&min_new, &max, sample, &params);

                    min_ratio = FLINT_MIN(min_ratio, min_old/min_new);
                    flint_printf(" %.2f ", min_old/min_new);
                    fflush(stdout);
                }

                flint_printf(") min %0.2f\n", min_ratio);

                /* assume that blas gets faster with more threads */
                if (min_ratio > 1)
                    break;
            }
        }
    }

    return 0;
}
