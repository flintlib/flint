/*
    Copyright 2009 William Hart
    Copyright 2010, 2021 Fredrik Johansson
    Copyright 2020 Daniel Schultz

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
#include "thread_support.h"
#include "perm.h"

typedef struct
{
    slong n;
    mp_limb_t modulus;
    int algorithm;
} mat_lu_t;

void sample(void * arg, ulong count)
{
    mat_lu_t * params = (mat_lu_t *) arg;
    int algorithm = params->algorithm;
    slong * P;
    nmod_mat_t A, LU;
    ulong i;
    flint_rand_t state;

    flint_randinit(state);

    nmod_mat_init(A, params->n, params->n, params->modulus);
    nmod_mat_init(LU, params->n, params->n, params->modulus);
    nmod_mat_randfull(A, state);
    P = _perm_init(params->n);

    prof_start();

    if (algorithm == 0)
        for (i = 0; i < count; i++)
        {
            nmod_mat_set(LU, A);
            nmod_mat_lu(P, LU, 0);
        }
    else if (algorithm == 1)
        for (i = 0; i < count; i++)
        {
            nmod_mat_set(LU, A);
            nmod_mat_lu_classical(P, LU, 0);
        }
    else if (algorithm == 2)
        for (i = 0; i < count; i++)
        {
            nmod_mat_set(LU, A);
            nmod_mat_lu_classical_delayed(P, LU, 0);
        }
    else
        for (i = 0; i < count; i++)
        {
            nmod_mat_set(LU, A);
            nmod_mat_lu_recursive(P, LU, 0);
        }

    prof_stop();

    nmod_mat_clear(A);
    nmod_mat_clear(LU);
    _perm_clear(P);

    flint_randclear(state);
}

slong bits_tab[] = { 5, 14, 15, 25, 30, 31, 32, 33, 60, 61, 62, 63, 64, 0 };

int main(void)
{
    double max;
    mat_lu_t params;
    slong i, bits, n;

    flint_printf("nmod_mat_lu:\n");

    flint_set_num_threads(1);
    flint_printf("threads = %wd\n", flint_get_num_threads());

    for (i = 0; (bits = bits_tab[i]) != 0; i++)
    {
        for (n = 4; n <= 1000; n += n/4 + 1)
        {
            double min_default = 0, min_classical = 0, min_delayed = 0, min_recursive = 0;

            params.n = n;
            params.modulus = n_nextprime(UWORD(1) << (bits - 1), 0);

            params.algorithm = 0;
            prof_repeat(&min_default, &max, sample, &params);

            params.algorithm = 1;
            prof_repeat(&min_classical, &max, sample, &params);

            params.algorithm = 2;
            prof_repeat(&min_delayed, &max, sample, &params);

            params.algorithm = 3;
            prof_repeat(&min_recursive, &max, sample, &params);

            flint_printf("b = %wd  n = %wd  default %.2f us  classical %.2f us  delayed %.2f us  recursive %.2f us\n",
                bits, n, min_default, min_classical, min_delayed, min_recursive);
        }
    }

    return 0;
}
