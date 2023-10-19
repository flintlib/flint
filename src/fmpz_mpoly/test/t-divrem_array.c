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

TEST_FUNCTION_START(fmpz_mpoly_divrem_array, state)
{
    int i, j, result, ok1, ok2;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r;
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
       fmpz_mpoly_init(r, ctx);

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
          fmpz_mpoly_randtest_bound(r, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok1 = fmpz_mpoly_divrem_array(k, r, h, g, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r, g, ctx);

          result = (ok1 == 0) || (ok1 && fmpz_mpoly_equal(f, k, ctx));

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
       fmpz_mpoly_clear(r, ctx);
    }

    /* Check f = g*q + r for random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r;
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
       fmpz_mpoly_init(r, ctx);

       len = n_randint(state, 16);
       len1 = n_randint(state, 16);
       len2 = n_randint(state, 16) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r, ctx);
          if (ok1)
		    {
             fmpz_mpoly_remainder_test(r, g, ctx);
             fmpz_mpoly_mul_johnson(k, h, g, ctx);
			    fmpz_mpoly_add(k, k, r, ctx);
          }

          result = (ok1 == 0) || (ok1 == 1 && fmpz_mpoly_equal(f, k, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check f = g*q + r for random polys\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
       fmpz_mpoly_clear(r, ctx);
    }

    /* Check aliasing of quotient with first argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, r1, r2;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(r1, ctx);
       fmpz_mpoly_init(r2, ctx);

       len = n_randint(state, 16);
       len1 = n_randint(state, 16);
       len2 = n_randint(state, 16) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(r2, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(f, r2, f, g, ctx);
          fmpz_mpoly_assert_canonical(f, ctx);
          fmpz_mpoly_assert_canonical(r2, ctx);
          if (ok2)
             fmpz_mpoly_remainder_test(r2, g, ctx);

          result = (ok1 == 0 || ok2 == 0) ||
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, f, ctx)
				                         && fmpz_mpoly_equal(r1, r2, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of quotient with first argument\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(r1, ctx);
       fmpz_mpoly_clear(r2, ctx);
    }

    /* Check aliasing of quotient with second argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, r1, r2;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(r1, ctx);
       fmpz_mpoly_init(r2, ctx);

       len = n_randint(state, 16);
       len1 = n_randint(state, 16);
       len2 = n_randint(state, 16) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(r2, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(g, r2, f, g, ctx);
          fmpz_mpoly_assert_canonical(g, ctx);
          fmpz_mpoly_assert_canonical(r2, ctx);

          result = (ok1 == 0 || ok2 == 0) ||
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, g, ctx)
				                         && fmpz_mpoly_equal(r1, r2, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of quotient with second argument\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(r1, ctx);
       fmpz_mpoly_clear(r2, ctx);
    }

    /* Check aliasing of remainder with first argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r1;
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
       fmpz_mpoly_init(r1, ctx);

       len = n_randint(state, 16);
       len1 = n_randint(state, 16);
       len2 = n_randint(state, 16) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(k, f, f, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);
          fmpz_mpoly_assert_canonical(f, ctx);
          if (ok2)
             fmpz_mpoly_remainder_test(f, g, ctx);

          result = (ok1 == 0 || ok2 == 0) ||
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, k, ctx)
				                         && fmpz_mpoly_equal(r1, f, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of remainder with first argument\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
       fmpz_mpoly_clear(r1, ctx);
    }

    /* Check aliasing of remainder with second argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r1;
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
       fmpz_mpoly_init(r1, ctx);

       len = n_randint(state, 16);
       len1 = n_randint(state, 16);
       len2 = n_randint(state, 16) + 1;

       exp_bound =  n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(k, g, f, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);
          fmpz_mpoly_assert_canonical(g, ctx);

          result = (ok1 == 0 || ok2 == 0) ||
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, k, ctx)
				                         && fmpz_mpoly_equal(r1, g, ctx));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of remainder with second argument\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
       fmpz_mpoly_clear(r1, ctx);
    }

    TEST_FUNCTION_END(state);
}
