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
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randlimb(state);  // a is arbitrary

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randint(state, d);  // must be < d

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
        {
            const ulong aj_pr = n_mulmod_precomp_shoup(array[j], d);
            array[j] = n_mulmod_shoup(array[j], a, aj_pr, d);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_no_precomp(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, d);  // a must be < d

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);  // array[j] is arbitrary

        const ulong a_pr = n_mulmod_precomp_shoup(a, d);

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_mulmod_shoup(a, array[j], a_pr, d);
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_precomp_only(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong d = n_randbits(state, bits);  // 0 < d < 2**(FLINT_BITS-1) required by mulmod_shoup

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randint(state, d);  // must be < d

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
        {
            const ulong FLINT_SET_BUT_UNUSED(aj_pr) = n_mulmod_precomp_shoup(array[j], d);
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    double min, max;

    flint_printf("mulmod_shoup min time / max time is:\n");

    prof_repeat(&min, &max, sample, NULL);
    flint_printf("   - including precomputation: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_no_precomp, NULL);
    flint_printf("   - excluding precomputation: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_precomp_only, NULL);
    flint_printf("   - precomputation alone: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    return 0;
}
