/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>

#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result, max_threads = 5;
    int tmul = 20;
    FLINT_TEST_INIT(state);
#ifdef _WIN32
    tmul = 1;
#endif

    flint_printf("mul_array_threaded....");
    fflush(stdout);

    /* Check mul_array_threaded matches mul_johnson */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, max_bound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        max_bound = ctx->minfo->ord == ORD_LEX ? 200 : 100;
        max_bound = max_bound/ctx->minfo->nvars/ctx->minfo->nvars;
        exp_bound =  n_randint(state, max_bound) + 1;
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            result = fmpz_mpoly_mul_array_threaded(k, f, g, ctx);
            if (!result)
            {
                continue;
            }
            fmpz_mpoly_assert_canonical(k, ctx);
            result = fmpz_mpoly_equal(h, k, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check mul_array_threaded matches mul_johnson\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, max_bound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        max_bound = 200/ctx->minfo->nvars/ctx->minfo->nvars;
        exp_bound =  n_randint(state, max_bound) + 1;
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            result = fmpz_mpoly_mul_array_threaded(f, f, g, ctx);
            if (!result)
                continue;

            fmpz_mpoly_assert_canonical(f, ctx);
            result = fmpz_mpoly_equal(h, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, max_bound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        max_bound = 200/ctx->minfo->nvars/ctx->minfo->nvars;
        exp_bound =  n_randint(state, max_bound) + 1;
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            result = fmpz_mpoly_mul_array_threaded(g, f, g, ctx);
            if (!result)
                continue;

            fmpz_mpoly_assert_canonical(g, ctx);
            result = fmpz_mpoly_equal(h, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second argument\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

