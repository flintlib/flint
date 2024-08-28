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

#define NB_ITER 1000  // multiple of 4

void sample(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    nn_ptr array_pr = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong n = n_randbits(state, bits);  // 0 < n < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, n);  // a < n required for precomputation
        ulong a_pr_quo, a_pr_rem;
        n_mulmod_precomp_shoup_quo_rem(&a_pr_quo, &a_pr_rem, a, n);

        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j] = n_randint(state, n);  // < n required for precomputation
            array_pr[j] = n_mulmod_precomp_shoup(array[j], n);
        }

        prof_start();
        // note: unrolling does not seem to help
        for (ulong j = 0; j < NB_ITER; j++)
        {
            n_mulmod_and_precomp_shoup(array+j+0, array_pr+j+0, a, array[j+0], a_pr_quo, a_pr_rem, array_pr[j+0], n);
            //n_mulmod_and_precomp_shoup(array+j+1, array_pr+j+1, a, array[j+1], a_pr_quo, a_pr_rem, array_pr[j+1], n);
            //n_mulmod_and_precomp_shoup(array+j+2, array_pr+j+2, a, array[j+2], a_pr_quo, a_pr_rem, array_pr[j+2], n);
            //n_mulmod_and_precomp_shoup(array+j+3, array_pr+j+3, a, array[j+3], a_pr_quo, a_pr_rem, array_pr[j+3], n);
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
        ulong n = n_randbits(state, bits);  // 0 < n < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, n);  // a < d required for precomputation
        ulong a_precomp = n_mulmod_precomp_shoup(a, n);

        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j] = n_randint(state, n);  // < n required for precomputation
        }

        prof_start();
        // note: unrolling does not seem to help
        for (ulong j = 0; j < NB_ITER; j++)
        {
            array[j+0] = n_mulmod_shoup(a, array[j+0], a_precomp, n);
            array_pr[j+0] = n_mulmod_precomp_shoup(array[j+0], n);
            //array[j+1] = n_mulmod_shoup(a, array[j+1], a_precomp, n);
            //array_pr[j+1] = n_mulmod_precomp_shoup(array[j+1], n);
            //array[j+2] = n_mulmod_shoup(a, array[j+2], a_precomp, n);
            //array_pr[j+2] = n_mulmod_precomp_shoup(array[j+2], n);
            //array[j+3] = n_mulmod_shoup(a, array[j+3], a_precomp, n);
            //array_pr[j+3] = n_mulmod_precomp_shoup(array[j+3], n);
        }
        prof_stop();
    }

    flint_free(array);
    flint_free(array_pr);
    FLINT_TEST_CLEAR(state);
}

/* given a in [0..n), compute the array [array[j], j=0...2*len-1] 
 * where array[2*j] = a**(j+1) mod n
 *   and array[2*j+1] = associated precomputed quotient
 **/
void sample_geometric(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(2*NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong n = n_randbits(state, bits);  // 0 < n < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, n);  // a < n required for precomputation

        prof_start();
        {
            ulong a_pr_quo, a_pr_rem;
            n_mulmod_precomp_shoup_quo_rem(&a_pr_quo, &a_pr_rem, a, n);
            // note: unrolling helps
            // note: computing with a**4 to reduce dependencies helps
            array[0] = a; array[1] = a_pr_quo;
            n_mulmod_and_precomp_shoup(array+(2*0+2), array+(2*0+3), a, array[2*0+0], a_pr_quo, a_pr_rem, array[2*0+1], n);
            n_mulmod_and_precomp_shoup(array+(2*0+4), array+(2*0+5), a, array[2*0+2], a_pr_quo, a_pr_rem, array[2*0+3], n);
            n_mulmod_and_precomp_shoup(array+(2*0+6), array+(2*0+7), a, array[2*0+4], a_pr_quo, a_pr_rem, array[2*0+5], n);
            // a**4 and its precomputed data
            a = array[6];
            a_pr_quo = array[7]; a_pr_rem = n_mulmod_precomp_shoup_rem_from_quo(a_pr_quo, n);
            for (ulong j = 4; j+3 < NB_ITER; j+=4)
            {
                n_mulmod_and_precomp_shoup(array+(2*j+0), array+(2*j+1), a, array[2*j-8], a_pr_quo, a_pr_rem, array[2*j-7], n);
                n_mulmod_and_precomp_shoup(array+(2*j+2), array+(2*j+3), a, array[2*j-6], a_pr_quo, a_pr_rem, array[2*j-5], n);
                n_mulmod_and_precomp_shoup(array+(2*j+4), array+(2*j+5), a, array[2*j-4], a_pr_quo, a_pr_rem, array[2*j-3], n);
                n_mulmod_and_precomp_shoup(array+(2*j+6), array+(2*j+7), a, array[2*j-2], a_pr_quo, a_pr_rem, array[2*j-1], n);
            }
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_geometric_direct(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(2*NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong n = n_randbits(state, bits);  // 0 < n < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, n);  // a < n required for precomputation

        prof_start();
        {
            const ulong a_precomp = n_mulmod_precomp_shoup(a, n);
            array[0] = a; array[1] = a_precomp;
            array[2*0+2] = n_mulmod_shoup(a, array[2*0+0], a_precomp, n);
            array[2*0+3] = n_mulmod_precomp_shoup(array[2*0+2], n);
            array[2*0+4] = n_mulmod_shoup(a, array[2*0+2], a_precomp, n);
            array[2*0+5] = n_mulmod_precomp_shoup(array[2*0+4], n);
            array[2*0+6] = n_mulmod_shoup(a, array[2*0+4], a_precomp, n);
            array[2*0+7] = n_mulmod_precomp_shoup(array[2*0+6], n);
            // note: unrolling helps a bit
            // note: computing with a**4 does not help
            for (ulong j = 4; j+3 < NB_ITER; j+=4)
            {
                array[2*j+0] = n_mulmod_shoup(a, array[2*j-2], a_precomp, n);
                array[2*j+1] = n_mulmod_precomp_shoup(array[2*j+0], n);
                array[2*j+2] = n_mulmod_shoup(a, array[2*j+0], a_precomp, n);
                array[2*j+3] = n_mulmod_precomp_shoup(array[2*j+2], n);
                array[2*j+4] = n_mulmod_shoup(a, array[2*j+2], a_precomp, n);
                array[2*j+5] = n_mulmod_precomp_shoup(array[2*j+4], n);
                array[2*j+6] = n_mulmod_shoup(a, array[2*j+4], a_precomp, n);
                array[2*j+7] = n_mulmod_precomp_shoup(array[2*j+6], n);
            }
        }
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

/* given a in [0..n), compute the arrays
 *  [a**(j+1) mod n,  j=0...] and [associated precomputed quotients]
 *  --> similar speed as single array variant
 **/
void sample_geometric_twoarrays(void * arg, ulong count)
{
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    nn_ptr array_pr = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong bits = n_randint(state, FLINT_BITS - 1) + 1;  // 1...63
        ulong n = n_randbits(state, bits);  // 0 < n < 2**(FLINT_BITS-1) required by mulmod_shoup
        ulong a = n_randint(state, n);  // a < n required for precomputation

        prof_start();
        {
            ulong a_pr_quo, a_pr_rem;
            n_mulmod_precomp_shoup_quo_rem(&a_pr_quo, &a_pr_rem, a, n);
            array[0] = a;
            array_pr[0] = a_pr_quo;
            n_mulmod_and_precomp_shoup(array+0+1, array_pr+0+1, a, array[0+0], a_pr_quo, a_pr_rem, array_pr[0+0], n);
            n_mulmod_and_precomp_shoup(array+0+2, array_pr+0+2, a, array[0+1], a_pr_quo, a_pr_rem, array_pr[0+1], n);
            n_mulmod_and_precomp_shoup(array+0+3, array_pr+0+3, a, array[0+2], a_pr_quo, a_pr_rem, array_pr[0+2], n);
            a = array[3]; a_pr_quo = array_pr[3]; a_pr_rem = n_mulmod_precomp_shoup_rem_from_quo(a_pr_quo, n);
            for (ulong j = 4; j+3 < NB_ITER; j+=4)
            {
                n_mulmod_and_precomp_shoup(array+j+0, array_pr+j+0, a, array[j-4], a_pr_quo, a_pr_rem, array_pr[j-4], n);
                n_mulmod_and_precomp_shoup(array+j+1, array_pr+j+1, a, array[j-3], a_pr_quo, a_pr_rem, array_pr[j-3], n);
                n_mulmod_and_precomp_shoup(array+j+2, array_pr+j+2, a, array[j-2], a_pr_quo, a_pr_rem, array_pr[j-2], n);
                n_mulmod_and_precomp_shoup(array+j+3, array_pr+j+3, a, array[j-1], a_pr_quo, a_pr_rem, array_pr[j-1], n);
            }
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
    flint_printf("   - mulmod_and_precomp_shoup (excluding precomputation for the input): %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_direct, NULL);
    flint_printf("   - mulmod_shoup+precomp_shoup (excluding precomputation for the input): %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_geometric, NULL);
    flint_printf("   - geometric with mulmod_and_precomp_shoup: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_geometric_direct, NULL);
    flint_printf("   - geometric with mulmod_shoup+precomp_shoup: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);

    prof_repeat(&min, &max, sample_geometric_twoarrays, NULL);
    flint_printf("   - geometric with mulmod_and_precomp_shoup, variant: %.3f cycles / %.3f cycles\n",
            (min/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER, (max/(double)FLINT_CLOCK_SCALE_FACTOR)/NB_ITER);


    return 0;
}
