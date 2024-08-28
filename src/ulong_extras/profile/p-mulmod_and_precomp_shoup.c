/*
   Copyright 2024 (C) Vincent Neiger

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
   */

#include "profiler.h"
#include "ulong_extras.h"

#define NB_ITER 1000

void sample(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    nn_ptr array_pr = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a < d required for precomputation
        ulong a_pr_quo, a_pr_rem;
        n_mulmod_precomp_shoup_quo_rem(&a_pr_quo, &a_pr_rem, a, d);

        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j] = n_randint(state, d);  // < d required for precomputation
            array_pr[j] = n_mulmod_precomp_shoup(array[j], d);
        }

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
        {
            n_mulmod_and_precomp_shoup(array+j, array_pr+j, a, array[j], a_pr_quo, a_pr_rem, array_pr[j], d);
        }
        prof_stop();
    }

    flint_free(array);
    flint_free(array_pr);
    FLINT_TEST_CLEAR(state);
}

void sample_direct(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    nn_ptr array_pr = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a < d required for precomputation
        ulong a_precomp = n_mulmod_precomp_shoup(a, d);

        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j] = n_randint(state, d);  // < d required for precomputation
        }

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j] = n_mulmod_shoup(a, array[j], a_precomp, d);
            array_pr[j] = n_mulmod_precomp_shoup(array[j], d);
        }
        prof_stop();
    }

    flint_free(array);
    flint_free(array_pr);
    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    double min, max;

    flint_printf("mulmod_shoup min time / max time is:\n");

    prof_repeat(&min, &max, sample, NULL);
    flint_printf("   - mulmod_and_precomp_shoup: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample, NULL);
    flint_printf("   - mulmod_shoup+precomp_shoup: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    return 0;
}
