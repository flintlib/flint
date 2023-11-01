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

TEST_FUNCTION_START(fmpz_mpoly_divides_array, state)
{
    int i, j, result, ok1, ok2;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok1 = fmpz_mpoly_divides_array(k, h, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);

          result = (ok1 == -1) || (ok1 && fmpz_mpoly_equal(f, k, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
    }

    /* Check random polys don't divide */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);

       len = n_randint(state, 30);
       len1 = n_randint(state, 30);
       len2 = n_randint(state, 30) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_divides_array(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          if (ok1)
          {
             fmpz_mpoly_mul_johnson(k, h, g, ctx);
             fmpz_mpoly_assert_canonical(k, ctx);
          }

          result = (ok1 == -1 || ok1 == 0) || (ok1 == 1 && fmpz_mpoly_equal(f, k, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check random polys don't divide\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok1 = fmpz_mpoly_divides_array(k, h, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);
          ok2 = fmpz_mpoly_divides_array(h, h, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          result = (ok1 == -1 || ok2 == -1) ||
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, k, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing first argument\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok1 = fmpz_mpoly_divides_array(k, h, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);
          ok2 = fmpz_mpoly_divides_array(g, h, g, ctx);
          fmpz_mpoly_assert_canonical(g, ctx);

          result = (ok1 == -1 && ok2 == -1) ||
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(g, k, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing second argument\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
    }

    TEST_FUNCTION_END(state);
}
