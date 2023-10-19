/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_mul, state)
{
    slong i, j, max_threads = 5;
    slong tmul = 10;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 2, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k1, ctx);
        fmpz_mod_mpoly_init(k2, ctx);
        fmpz_mod_mpoly_init(t1, ctx);
        fmpz_mod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        fmpz_mod_mpoly_randtest_bits(k1, state, len, exp_bits, ctx);
        fmpz_mod_mpoly_randtest_bits(k2, state, len, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fmpz_mod_mpoly_randtest_bits(h, state, len2, exp_bits, ctx);

            fmpz_mod_mpoly_add(t1, g, h, ctx);
            fmpz_mod_mpoly_assert_canonical(t1, ctx);
            fmpz_mod_mpoly_mul(k1, f, t1, ctx);
            fmpz_mod_mpoly_assert_canonical(k1, ctx);
            fmpz_mod_mpoly_mul(t1, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(t1, ctx);
            fmpz_mod_mpoly_mul(t2, f, h, ctx);
            fmpz_mod_mpoly_assert_canonical(t2, ctx);
            fmpz_mod_mpoly_add(k2, t1, t2, ctx);
            fmpz_mod_mpoly_assert_canonical(k2, ctx);

            if (!fmpz_mod_mpoly_equal(k1, k2, ctx))
            {
                flint_printf("FAIL: Check f*(g + h) = f*g + f*h\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k1, ctx);
        fmpz_mod_mpoly_clear(k2, ctx);
        fmpz_mod_mpoly_clear(t1, ctx);
        fmpz_mod_mpoly_clear(t2, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bound, exp_bound1, exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 4, 200);
        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);

        exp_bound = 3 + n_randint(state, 1 + 100/n/n);
        exp_bound1 = 3 + n_randint(state, 1 + 100/n/n);
        exp_bound2 = 3 + n_randint(state, 1 + 100/n/n);

        len = exp_bound + 1;
        len1 = exp_bound1 + 1;
        len2 = exp_bound2 + 1;
        for (j = n_randint(state, n) + 2; j >= 0; j--)
        {
            len *= exp_bound + 1;
            len1 *= exp_bound1 + 1;
            len2 *= exp_bound2 + 1;
            len = FLINT_MIN(len, WORD(100));
            len1 = FLINT_MIN(len, WORD(100));
            len2 = FLINT_MIN(len, WORD(100));
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            fmpz_mod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mod_mpoly_mul(h, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            fmpz_mod_mpoly_mul(f, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            if (!fmpz_mod_mpoly_equal(h, f, ctx))
            {
                flint_printf("FAIL: Check aliasing first arg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            fmpz_mod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mod_mpoly_mul(h, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            fmpz_mod_mpoly_mul(g, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(g, ctx);
            if (!fmpz_mod_mpoly_equal(h, g, ctx))
            {
                flint_printf("FAIL: Check aliasing second arg\n");
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
