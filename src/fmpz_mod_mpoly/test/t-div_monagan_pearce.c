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

TEST_FUNCTION_START(fmpz_mod_mpoly_div_monagan_pearce, state)
{
    slong i, j;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k, l;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 20, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);
        fmpz_mod_mpoly_init(l, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(l, state, len, exp_bits, ctx);

            fmpz_mod_mpoly_mul(h, f, g, ctx);

            fmpz_mod_mpoly_div_monagan_pearce(k, h, g, ctx);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            if (!fmpz_mod_mpoly_equal(k, f, ctx))
            {
                flint_printf("FAIL: Check f*g/g = f\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_set(l, h, ctx);
            fmpz_mod_mpoly_div_monagan_pearce(l, l, g, ctx);
            fmpz_mod_mpoly_assert_canonical(l, ctx);
            if (!fmpz_mod_mpoly_equal(l, f, ctx))
            {
                flint_printf("FAIL: Check f*g/g = f aliasing dividend\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_div_monagan_pearce(g, h, g, ctx);
            fmpz_mod_mpoly_assert_canonical(g, ctx);
            if (!fmpz_mod_mpoly_equal(g, f, ctx))
            {
                flint_printf("FAIL: Check f*g/g = f aliasing divisor\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_clear(l, ctx);

        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check div matches divrem for random polys */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, q, r, k;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 20, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(q, ctx);
        fmpz_mod_mpoly_init(r, ctx);
        fmpz_mod_mpoly_init(k, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 50/n) + 1;
        exp_bound1 = n_randint(state, 50/n) + 1;
        exp_bound2 = n_randint(state, 50/n) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bound(q, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_divrem_monagan_pearce(q, r, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(q, ctx);
            fmpz_mod_mpoly_assert_canonical(r, ctx);
            fmpz_mod_mpoly_remainder_strongtest(r, g, ctx);

            fmpz_mod_mpoly_mul(k, q, g, ctx);
            fmpz_mod_mpoly_add(k, k, r, ctx);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            if (!fmpz_mod_mpoly_equal(f, k, ctx))
            {
                flint_printf("FAIL: Check f = g*q + r for random polys\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_div_monagan_pearce(k, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            if (!fmpz_mod_mpoly_equal(k, q, ctx))
            {
                flint_printf("FAIL: Check div matches divrem\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_set(k, f, ctx);
            fmpz_mod_mpoly_div_monagan_pearce(k, k, g, ctx);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            if (!fmpz_mod_mpoly_equal(k, q, ctx))
            {
                flint_printf("FAIL: Check div matches divrem aliasing dividend\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_div_monagan_pearce(g, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(g, ctx);
            if (!fmpz_mod_mpoly_equal(g, q, ctx))
            {
                flint_printf("FAIL: Check div matches divrem aliasing divisor\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(q, ctx);
        fmpz_mod_mpoly_clear(r, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
