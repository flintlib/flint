/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("get_set_string_pretty....");
    fflush(stdout);

    {
        slong len1;
        flint_bitcnt_t exp_bits1;
        mp_limb_t modulus;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, f1;
        char * str;
        const char * vars[] = {"x","y","z","w","u","v"};

        for (i = 0; i < flint_test_multiplier(); i++)
        {
            modulus = n_randint(state, FLINT_BITS - 1) + 1;
            modulus = n_randbits(state, modulus);
            modulus = n_nextprime(modulus, 1);
            nmod_mpoly_ctx_init_rand(ctx, state, 6, modulus);

            nmod_mpoly_init(f, ctx);
            nmod_mpoly_init(f1, ctx);

            for (len1 = 3; len1 < 1000; len1 += len1/2)
            {
                exp_bits1 = 200;
                nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);

                str = nmod_mpoly_get_str_pretty(f, vars, ctx);
                nmod_mpoly_set_str_pretty(f1, str, vars, ctx);
                flint_free(str);

                if (!nmod_mpoly_equal(f, f1, ctx))
                {
                    flint_printf("FAIL\n");
                    fflush(stdout);
                    flint_abort();
                }

            }

            nmod_mpoly_clear(f, ctx);
            nmod_mpoly_clear(f1, ctx);

            nmod_mpoly_ctx_clear(ctx);
        }
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

