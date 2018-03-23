/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_mul_ui....");
    fflush(stdout);

    /* Check (f*a)*b = f*(a*b) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        ordering_t ord;
        mp_limb_t modulus;
        ulong a, b, c;
        slong nvars, len, exp_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init(ctx, nvars, ord, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            a = n_randbits(state, n_randint(state, FLINT_BITS/2) + 1);
            b = n_randbits(state, n_randint(state, FLINT_BITS/2) + 1);
            c = a*b;

            nmod_mpoly_scalar_mul_ui(g, f, a, ctx);
            nmod_mpoly_scalar_mul_ui(h, g, b, ctx);

            nmod_mpoly_scalar_mul_ui(k, f, c, ctx);

            result = nmod_mpoly_equal(h, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        ordering_t ord;
        mp_limb_t modulus;
        ulong c;
        slong nvars, len, exp_bits;

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

        exp_bits = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            c = n_randtest(state);

            nmod_mpoly_set(g, f, ctx);

            nmod_mpoly_scalar_mul_ui(h, f, c, ctx);

            nmod_mpoly_scalar_mul_ui(g, g, c, ctx);

            result = nmod_mpoly_equal(g, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

