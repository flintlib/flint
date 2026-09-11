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
#include "double_extras.h"

#define NB_ITER 1000

typedef struct
{
    ulong bits;
    ulong exp;
} info_t;


void sample_preinv(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;
    ulong bits = info->bits;
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong n = n_randbits(state, bits);  // 0 < n < 2**(FLINT_BITS)
        ulong ninv = n_preinvert_limb(n);
        ulong norm = flint_clz(n);

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randint(state, n);  // 0 <= array[j] < n

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_powmod_ui_preinv(array[j], exp, n, ninv, norm);
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_preinv2(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;
    ulong bits = info->bits;
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong n = n_randbits(state, bits);  // 0 < n < 2**(FLINT_BITS)
        ulong ninv = n_preinvert_limb(n);

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randlimb(state);

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_powmod2_ui_preinv(array[j], exp, n, ninv);
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

void sample_precomp(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;
    ulong bits = info->bits;
    nn_ptr array = (nn_ptr) flint_malloc(NB_ITER*sizeof(ulong));
    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        ulong n = n_randbits(state, bits);  // 0 < n < 2**bits
        double ninv = n_precompute_inverse(n);

        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_randint(state, n);  // 0 <= array[j] < n

        prof_start();
        for (ulong j = 0; j < NB_ITER; j++)
            array[j] = n_powmod_ui_precomp(array[j], exp, n, ninv);
        prof_stop();
    }

    flint_free(array);
    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    double min, max;

    const ulong bits_nb = 5;
    ulong bits_list[] = {20, 30, 50, 60, 64};
    const ulong exp_nb = 11;
    ulong exp_list[] = {5, 10, 20, 40, 80, 160, 1000, 10000, 100000, 1000000L, 10000000L};

    flint_printf("compute an exponentiation a**e mod n, with nbits(n) = b\n");
    flint_printf("  computation is repeated on the element of a %wu-length array\n");
    flint_printf("  time is divided by %wu * FLINT_CLOCK_SCALE_FACTOR * log_2(exp)\n", NB_ITER, NB_ITER);
    flint_printf("timings are: powmod_ui_precomp | powmod_ui_preinv | powmod2_ui_preinv\n");
    flint_printf("b \\ e\t");
    for (ulong e = 0; e < exp_nb; e++)
        flint_printf("%wu\t\t", exp_list[e]);
    flint_printf("\n");

    info_t info;

    for (ulong b = 0; b < bits_nb; b++)
    {
        info.bits = bits_list[b];
        flint_printf("%wu\t", info.bits);

        for (ulong e = 0; e < exp_nb; e++)
        {
            info.exp = exp_list[e];
            double log_exp = d_log2((double)info.exp);

            if (info.bits <= 53)
            {
                prof_repeat(&min, &max, sample_precomp, (void *) &info);
                flint_printf("%4.1f|", min/(NB_ITER * FLINT_CLOCK_SCALE_FACTOR * log_exp));
            }
            else
                flint_printf(" na |");

            prof_repeat(&min, &max, sample_preinv, (void *) &info);
            flint_printf("%4.1f|", min/(NB_ITER * FLINT_CLOCK_SCALE_FACTOR * log_exp));

            prof_repeat(&min, &max, sample_preinv2, (void *) &info);
            flint_printf("%4.1f\t", min/(NB_ITER * FLINT_CLOCK_SCALE_FACTOR * log_exp));
        }
        flint_printf("\n");
    }

    return 0;
}
