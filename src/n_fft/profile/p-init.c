/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod.h"
#include "profiler.h"
#include "n_fft.h"

#define num_primes 5

typedef struct
{
   ulong prime;
   ulong depth;
   ulong maxdepth;
} info_t;

void sample_init2_root(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    ulong p = info->prime;
    ulong depth = info->depth;
    ulong maxdepth = info->maxdepth;

    const ulong len = UWORD(1) << depth;
    const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

    // modulus, roots of unity
    nmod_t mod;
    nmod_init(&mod, p);
    ulong cofactor = (p - 1) >> maxdepth;
    ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), cofactor, mod);
    ulong w = nmod_pow_ui(w0, 1UL<<(maxdepth - depth), mod);

    FLINT_TEST_INIT(state);

    for (ulong i = 0; i < count; i++)
    {
        prof_start();
        for (ulong j = 0; j < rep; j++)
        {
            n_fft_ctx_t F;
            n_fft_ctx_init2_root(F, w, depth, cofactor, depth, p);
            n_fft_ctx_clear(F);
        }
        prof_stop();
    }

    FLINT_TEST_CLEAR(state);
}

/*-----------------------------------------------------------------*/
/* initialize context for FFT for several bit lengths and depths   */
/*-----------------------------------------------------------------*/
void time_fft_init(ulong * primes, ulong * max_depths)
{
    for (ulong depth = 3; depth <= 25; depth++)
    {
        printf("%ld\t", depth);
        for (ulong k = 0; k < num_primes; k++)
        {
            if (depth <= max_depths[k])
            {
                info_t info;
                info.prime = primes[k];
                info.maxdepth = max_depths[k];
                info.depth = depth;

                const ulong len = UWORD(1) << depth;
                const ulong rep = FLINT_MAX(1, FLINT_MIN(1000, 1000000/len));

                double min;
                double max;

                prof_repeat(&min, &max, sample_init2_root, (void *) &info);

                flint_printf("%.1e|%.1e\t",
                        min/(double)1000000/rep,
                        min/(double)FLINT_CLOCK_SCALE_FACTOR/len/rep
                        );
            }
            else
                flint_printf("  na   |  na   \t");
        }
        flint_printf("\n");
    }

}

/*------------------------------------------------------------*/
/* main just calls time_init_set()                            */
/*------------------------------------------------------------*/
int main()
{
    printf("- depth == precomputing w**k, 0 <= k < 2**depth\n");
    printf("- timing init FFT context + clear at this depth:\n");
    printf("      t_raw == raw time\n");
    printf("      t_unit == raw time divided by 2**depth * clock scale factor\n");
    printf("\n");

    printf("     \t    20 bits    \t    31 bits    \t    42 bits    \t    50 bits    \t    60 bits    \n");
    printf("depth\tt_raw  | t_unit\tt_raw  | t_unit\tt_raw  | t_unit\tt_raw  | t_unit\tt_raw  | t_unit\n");

    // TODO fix for FLINT_BITS==32
    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_depths[num_primes] = { 18, 25, 25, 25, 25 };

    time_fft_init(primes, max_depths);

    return 0;
}

