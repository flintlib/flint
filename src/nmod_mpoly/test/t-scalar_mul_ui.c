/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_scalar_mul_ui, state)
{
    int i, j, result;

    /* Check (f*a)*b = f*(a*b) */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        ulong a, b, c;
        slong len1, len2, len3, len4;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3, exp_bits4;
        mp_limb_t modulus;

        modulus = n_randbits(state, 1 + n_randint(state, FLINT_BITS));
        modulus = FLINT_MAX(UWORD(2), modulus);

        nmod_mpoly_ctx_init_rand(ctx, state, 2, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        len3 = n_randint(state, 100);
        len4 = n_randint(state, 100);
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;
        exp_bits3 = n_randint(state, 200) + 2;
        exp_bits4 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len3, exp_bits3, ctx);
            nmod_mpoly_randtest_bits(k, state, len4, exp_bits4, ctx);

            a = n_randtest(state);
            b = n_randtest(state);

            nmod_mpoly_scalar_mul_ui(g, f, a, ctx);
            nmod_mpoly_scalar_mul_ui(h, g, b, ctx);
            NMOD_RED(a, a, ctx->mod);
            NMOD_RED(b, b, ctx->mod);
            c = nmod_mul(a, b, ctx->mod);
            nmod_mpoly_scalar_mul_ui(k, f, c, ctx);
            result = nmod_mpoly_equal(h, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check (f*a)*b = f*(a*b)\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            nmod_mpoly_scalar_mul_ui(g, f, a, ctx);
            nmod_mpoly_scalar_mul_ui(g, g, b, ctx);
            nmod_mpoly_scalar_mul_ui(f, f, c, ctx);
            result = nmod_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check (f*a)*b = f*(a*b) with aliasing\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*a*inv(a) = f */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        ulong a, b;
        slong len1, len2, len3, len4;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3, exp_bits4;
        mp_limb_t modulus;

        modulus = n_randbits(state, 1 + n_randint(state, FLINT_BITS));
        modulus = FLINT_MAX(UWORD(2), modulus);

        nmod_mpoly_ctx_init_rand(ctx, state, 2, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        len3 = n_randint(state, 100);
        len4 = n_randint(state, 100);
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;
        exp_bits3 = n_randint(state, 200) + 2;
        exp_bits4 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len3, exp_bits3, ctx);
            nmod_mpoly_randtest_bits(k, state, len4, exp_bits4, ctx);

            a = n_randtest(state);

            nmod_mpoly_scalar_mul_ui(g, f, a, ctx);

            NMOD_RED(a, a, ctx->mod);
            if (n_gcd(a, ctx->mod.n) != UWORD(1))
                continue;

            b = nmod_inv(a, ctx->mod);
            nmod_mpoly_scalar_mul_ui(h, g, b, ctx);
            result = nmod_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*a*inv(a) = f\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            nmod_mpoly_scalar_mul_ui(g, f, a, ctx);
            nmod_mpoly_scalar_mul_ui(g, g, b, ctx);
            result = nmod_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*a*inv(a) = f with aliasing\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
