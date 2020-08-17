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
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("add_sub_ui....");
    fflush(stdout);

    /* Check (f + c) - c = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        ordering_t ord;
        mp_limb_t modulus;
        ulong c;
        slong nvars, len, len1, len2;
        slong exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            c = n_randtest(state);

            nmod_mpoly_add_ui(g, f, c, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            nmod_mpoly_sub_ui(h, g, c, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = nmod_mpoly_equal(f, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check (f + c) - c = f\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g;
        ordering_t ord;
        mp_limb_t modulus;
        ulong c;
        slong nvars, len1, len2;
        slong exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_set(g, f, ctx);
            c = n_randtest(state);
            nmod_mpoly_add_ui(f, f, c, ctx);
            nmod_mpoly_assert_canonical(f, ctx);
            nmod_mpoly_sub_ui(f, f, c, ctx);
            nmod_mpoly_assert_canonical(f, ctx);
            result = nmod_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }



    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

