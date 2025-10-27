/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

/* helper to create examples */
void _nmod_vec_rand_not_zero(nn_ptr vec, flint_rand_t state, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
        vec[k] = 1 + n_randint(state, mod.n - 1);
}

typedef struct
{
   ulong prime;
   slong length;
} info_t;

void sample_invert_naive(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    nmod_t mod;
    nmod_init(&mod, info->prime);
    const ulong len = info->length;

    nn_ptr vec = (nn_ptr) flint_malloc(len*sizeof(ulong));
    nn_ptr res = _nmod_vec_init(len);
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        _nmod_vec_rand_not_zero(vec, state, len, mod);

        prof_start();
        _nmod_vec_invert_naive(res, vec, len, mod);
        prof_stop();
    }

    flint_free(vec);
    flint_free(res);
    FLINT_TEST_CLEAR(state);
}

void sample_invert_generic(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    nmod_t mod;
    nmod_init(&mod, info->prime);
    const ulong len = info->length;

    nn_ptr vec = (nn_ptr) flint_malloc(len*sizeof(ulong));
    nn_ptr res = _nmod_vec_init(len);
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        _nmod_vec_rand_not_zero(vec, state, len, mod);

        prof_start();
        _nmod_vec_invert_generic(res, vec, len, mod);
        prof_stop();
    }

    flint_free(vec);
    flint_free(res);
    FLINT_TEST_CLEAR(state);
}

void sample_invert_shoup(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    nmod_t mod;
    nmod_init(&mod, info->prime);
    const ulong len = info->length;

    nn_ptr vec = (nn_ptr) flint_malloc(len*sizeof(ulong));
    nn_ptr res = _nmod_vec_init(len);
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        _nmod_vec_rand_not_zero(vec, state, len, mod);

        prof_start();
        _nmod_vec_invert_shoup(res, vec, len, mod);
        prof_stop();
    }

    flint_free(vec);
    flint_free(res);
    FLINT_TEST_CLEAR(state);
}


int main(void)
{
    FLINT_TEST_INIT(state);

    double min, max;
    double mins_naive[18]; // note: max seems to be consistently identical or extremely close to min
    double mins_shoup[18];
    double mins_generic[18];
    info_t info;
    flint_bitcnt_t i;

    flint_printf("unit: all measurements in c/l (up to constant multiplicative factor)\n");
    flint_printf("profiled: naive | precomp shoup | generic\n");
    flint_printf("bit/len\t");
    for (int len = 1; len <= 16; ++len)
        flint_printf("%d\t\t", len);
    flint_printf("1024\t\t");
    flint_printf("65536\n");

    for (i = 2; i <= FLINT_BITS; i++)
    {
        info.prime = n_randprime(state, i, 1);

        for (int len = 1; len <= 16; ++len)
        {
            info.length = len;

            prof_repeat(&min, &max, sample_invert_naive, (void *) &info);
            mins_naive[len-1] = min;
            prof_repeat(&min, &max, sample_invert_shoup, (void *) &info);
            mins_shoup[len-1] = min;
            prof_repeat(&min, &max, sample_invert_generic, (void *) &info);
            mins_generic[len-1] = min;
        }

        info.length = 1024;
        prof_repeat(&min, &max, sample_invert_naive, (void *) &info);
        mins_naive[16] = min;
        prof_repeat(&min, &max, sample_invert_shoup, (void *) &info);
        mins_shoup[16] = min;
        prof_repeat(&min, &max, sample_invert_generic, (void *) &info);
        mins_generic[16] = min;

        info.length = 65536;
        prof_repeat(&min, &max, sample_invert_naive, (void *) &info);
        mins_naive[17] = min;
        prof_repeat(&min, &max, sample_invert_shoup, (void *) &info);
        mins_shoup[17] = min;
        prof_repeat(&min, &max, sample_invert_generic, (void *) &info);
        mins_generic[17] = min;

        if (i < FLINT_BITS)
        {
            flint_printf("%wd", i);
            for (int len = 1; len <= 16; ++len)
                flint_printf("\t%.2lf|%.2lf|%.2lf",
                             (mins_naive[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100),
                             (mins_shoup[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100),
                             (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100));
            flint_printf("\t%.2lf|%.2lf|%.2lf",
                         (mins_naive[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100),
                         (mins_shoup[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100),
                         (mins_generic[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100));
            flint_printf("\t%.2lf|%.2lf|%.2lf",
                         (mins_naive[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100),
                         (mins_shoup[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100),
                         (mins_generic[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100));
            flint_printf("\n");
        }
        else
        {
            flint_printf("%wd", i);
            for (int len = 1; len <= 16; ++len)
                flint_printf("\t%.2lf| na |%.2lf",
                             (mins_naive[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100),
                             (mins_generic[len-1]/(double)FLINT_CLOCK_SCALE_FACTOR)/(len*100));
            flint_printf("\t%.2lf| na |%.2lf",
                         (mins_naive[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100),
                         (mins_generic[16]/(double)FLINT_CLOCK_SCALE_FACTOR)/(1024*100));
            flint_printf("\t%.2lf| na |%.2lf",
                         (mins_naive[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100),
                         (mins_generic[17]/(double)FLINT_CLOCK_SCALE_FACTOR)/(65536*100));
            flint_printf("\n");
        }
    }

    FLINT_TEST_CLEAR(state);
    return 0;
}
