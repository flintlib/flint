/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_mul_johnson, state)
{
    int i, j, result;

    /* Check mul_johnson matches mul_array */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
       slong n, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;
       int ok1;

       fmpz_mpoly_ctx_init_rand(ctx, state, 20);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
       exp_bound =  n_randint(state, 800/n/n) + 1;
       exp_bound1 = n_randint(state, 800/n/n) + 1;
       exp_bound2 = n_randint(state, 800/n/n) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_mul_array(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          fmpz_mpoly_mul_johnson(k, f, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);

          result = ok1 == 0 || fmpz_mpoly_equal(h, k, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check mul_johnson matches mul_array\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
    }

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest_bits(k1, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(k2, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_add(t1, g, h, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul_johnson(k1, f, t1, ctx);
            fmpz_mpoly_assert_canonical(k1, ctx);
            fmpz_mpoly_mul_johnson(t1, f, g, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul_johnson(t2, f, h, ctx);
            fmpz_mpoly_assert_canonical(t2, ctx);
            fmpz_mpoly_add(k2, t1, t2, ctx);
            fmpz_mpoly_assert_canonical(k2, ctx);
            result = fmpz_mpoly_equal(k1, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k1, ctx);
        fmpz_mpoly_clear(k2, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_mul_johnson(f, f, g, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            result = fmpz_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_mul_johnson(g, f, g, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            result = fmpz_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
