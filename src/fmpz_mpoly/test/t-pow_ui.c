/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

/* Defined in t-pow_fps.c and t-pow_ui.c */
#ifndef fmpz_mpoly_pow_naive
#define fmpz_mpoly_pow_naive fmpz_mpoly_pow_naive
void fmpz_mpoly_pow_naive(fmpz_mpoly_t res, fmpz_mpoly_t f,
                                                 slong n, fmpz_mpoly_ctx_t ctx)
{
   if (n == 0)
      fmpz_mpoly_one(res, ctx);
   else if (f->length == 0)
      fmpz_mpoly_zero(res, ctx);
   else if (n == 1)
      fmpz_mpoly_set(res, f, ctx);
   else
   {
      slong i;
      fmpz_mpoly_t pow;

      fmpz_mpoly_init(pow, ctx);
      fmpz_mpoly_set(pow, f, ctx);

      for (i = 1; i < n - 1; i++)
         fmpz_mpoly_mul_johnson(pow, pow, f, ctx);

      fmpz_mpoly_mul_johnson(res, pow, f, ctx);

      fmpz_mpoly_clear(pow, ctx);
   }
}
#endif

TEST_FUNCTION_START(fmpz_mpoly_pow_ui, state)
{
    slong i, j, tmul = 20;

    /* Check pow_ui against pow_naive */
    for (i = 0; i < 10 * tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1;
        ulong pow;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);

        pow = n_randint(state, 1 + 50/(len1 + 2));

        exp_bits = n_randint(state, 600) + 2;
        exp_bits1 = n_randint(state, 600) + 10;
        exp_bits1 = n_randint(state, exp_bits1) + 2; /* increase chances of lower values */

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_pow_ui(g, f, pow, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            fmpz_mpoly_pow_naive(h, f, pow, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            if (!fmpz_mpoly_equal(g, h, ctx))
            {
                flint_printf("FAIL: Check pow_ui against pow_naive\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        slong len, len1;
        ulong pow;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);

        pow = n_randint(state, 8);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);

        exp_bits = n_randint(state, 600) + 2;
        exp_bits1 = n_randint(state, 600) + 10;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_pow_ui(g, f, pow, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            fmpz_mpoly_pow_ui(f, f, pow, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            if (!fmpz_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL: Check aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
