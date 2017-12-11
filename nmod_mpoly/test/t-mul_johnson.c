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

    flint_printf("mul_johnson....");
    fflush(stdout);

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, t1, t2, k1, k2;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
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
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);
        nmod_mpoly_init(k1, ctx);
        nmod_mpoly_init(k2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        maxbits = FLINT_BITS - mpoly_ordering_isdeg(ctx->minfo)*FLINT_BIT_COUNT(nvars);
        exp_bits = n_randint(state, maxbits - 2) + 1;
        exp_bits1 = n_randint(state, maxbits - 2) + 1;
        exp_bits2 = n_randint(state, maxbits - 2) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);

            nmod_mpoly_add(t1, g, h, ctx);
            nmod_mpoly_test(t1, ctx);
            nmod_mpoly_mul_johnson(k1, f, t1, ctx);
            nmod_mpoly_test(k1, ctx);

            nmod_mpoly_mul_johnson(t1, f, g, ctx);
            nmod_mpoly_test(t1, ctx);
            nmod_mpoly_mul_johnson(t2, f, h, ctx);
            nmod_mpoly_test(t2, ctx);
            nmod_mpoly_add(k2, t1, t2, ctx);
            nmod_mpoly_test(k2, ctx);
            result = nmod_mpoly_equal(k1, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);
        nmod_mpoly_clear(k1, ctx);
        nmod_mpoly_clear(k2, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        ordering_t ord;
        mp_limb_t modulus;
        slong maxbits;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
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

        maxbits = FLINT_BITS - mpoly_ordering_isdeg(ctx->minfo)*FLINT_BIT_COUNT(nvars);
        exp_bits = n_randint(state, maxbits - 1) + 1;
        exp_bits1 = n_randint(state, maxbits - 2) + 1;
        exp_bits2 = n_randint(state, maxbits - 2) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);
            nmod_mpoly_mul_johnson(f, f, g, ctx);
            nmod_mpoly_test(f, ctx);
            result = nmod_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        ordering_t ord;
        mp_limb_t modulus;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
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

        exp_bits = n_randint(state, FLINT_BITS -
                     mpoly_ordering_isdeg(ctx->minfo)*FLINT_BIT_COUNT(nvars) - 1) + 1;
        exp_bits1 = n_randint(state, FLINT_BITS -
                     mpoly_ordering_isdeg(ctx->minfo)*FLINT_BIT_COUNT(nvars) - 1) + 1;
        exp_bits2 = n_randint(state, FLINT_BITS -
                     mpoly_ordering_isdeg(ctx->minfo)*FLINT_BIT_COUNT(nvars) - 1) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest(h, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_test(h, ctx);
            nmod_mpoly_mul_johnson(f, g, f, ctx);
            nmod_mpoly_test(f, ctx);
            result = nmod_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument\ni = %wd, j = %wd\n", i ,j);
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

