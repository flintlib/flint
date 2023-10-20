/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_scalar_mul_fq_nmod, state)
{
    int i, j, result;

    /* Check (f*a)*b = f*(a*b) */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k;
        fq_nmod_t a, b, c;
        slong len1, len2, len3, len4;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3, exp_bits4;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_init(a, ctx->fqctx);
        fq_nmod_init(b, ctx->fqctx);
        fq_nmod_init(c, ctx->fqctx);

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
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len3, exp_bits3, ctx);
            fq_nmod_mpoly_randtest_bits(k, state, len4, exp_bits4, ctx);

            fq_nmod_randtest(a, state, ctx->fqctx);
            fq_nmod_randtest(b, state, ctx->fqctx);
            fq_nmod_mul(c, a, b, ctx->fqctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(g, f, a, ctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(h, g, b, ctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(k, f, c, ctx);
            result = fq_nmod_mpoly_equal(h, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check (f*a)*b = f*(a*b)\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            fq_nmod_mpoly_scalar_mul_fq_nmod(g, f, a, ctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(g, g, b, ctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(f, f, c, ctx);
            result = fq_nmod_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check (f*a)*b = f*(a*b) with aliasing\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(k, ctx);
        fq_nmod_clear(a, ctx->fqctx);
        fq_nmod_clear(b, ctx->fqctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*a*inv(a) = f */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k;
        fq_nmod_t a, b;
        slong len1, len2, len3, len4;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3, exp_bits4;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_init(a, ctx->fqctx);
        fq_nmod_init(b, ctx->fqctx);

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
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len3, exp_bits3, ctx);
            fq_nmod_mpoly_randtest_bits(k, state, len4, exp_bits4, ctx);

            fq_nmod_randtest_not_zero(a, state, ctx->fqctx);
            fq_nmod_inv(b, a, ctx->fqctx);

            fq_nmod_mpoly_scalar_mul_fq_nmod(g, f, a, ctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(h, g, b, ctx);
            result = fq_nmod_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*a*inv(a) = f\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            fq_nmod_mpoly_scalar_mul_fq_nmod(g, f, a, ctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(g, g, b, ctx);
            result = fq_nmod_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*a*inv(a) = f with aliasing\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(k, ctx);
        fq_nmod_clear(a, ctx->fqctx);
        fq_nmod_clear(b, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
