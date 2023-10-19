/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_mul_dense, state)
{
    slong i, j;

    /* Check mul_dense matches mul_johnson */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        slong max_bound, exp_bound, exp_bound1, exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);

        len = n_randint(state, 200);
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 100/n/n;
        exp_bound = UWORD(1) << (FLINT_BITS - 1);
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            fmpz_mod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            if (!fmpz_mod_mpoly_mul_dense(k, f, g, ctx))
                continue;
            fmpz_mod_mpoly_assert_canonical(k, ctx);

            if (!fmpz_mod_mpoly_equal(h, k, ctx))
            {
                flint_printf("FAIL: Check mul_dense matches mul_johnson\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h;
        slong len, len1, len2;
        slong max_bound, exp_bound, exp_bound1, exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);

        len = n_randint(state, 200);
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 100/n/n;
        exp_bound = UWORD(1) << (FLINT_BITS - 1);
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            fmpz_mod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            if (!fmpz_mod_mpoly_mul_dense(f, f, g, ctx))
                continue;
            fmpz_mod_mpoly_assert_canonical(f, ctx);

            if (!fmpz_mod_mpoly_equal(h, f, ctx))
            {
                flint_printf("FAIL: Check aliasing first argument\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h;
        slong len, len1, len2;
        slong max_bound, exp_bound, exp_bound1, exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);

        len = n_randint(state, 200);
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 100/n/n;
        exp_bound = UWORD(1) << (FLINT_BITS - 1);
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            fmpz_mod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            if (!fmpz_mod_mpoly_mul_dense(f, g, f, ctx))
                continue;
            fmpz_mod_mpoly_assert_canonical(f, ctx);

            if (!fmpz_mod_mpoly_equal(h, f, ctx))
            {
                flint_printf("FAIL: Check aliasing first argument\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
