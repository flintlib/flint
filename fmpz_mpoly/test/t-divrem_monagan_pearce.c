/*
    Copyright (C) 2017 William Hart
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

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("divrem_monagan_pearce....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r;
       ordering_t ord;
       slong nvars, len, len1, len2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

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

       exp_bits = n_randint(state, 200) + 1;
       exp_bits1 = n_randint(state, 200) + 1;
       exp_bits2 = n_randint(state, 200) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
          do {
             fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
          fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);
          fmpz_mpoly_randtest_bits(r, state, len, coeff_bits, exp_bits, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);

          fmpz_mpoly_divrem_monagan_pearce(k, r, h, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);
          fmpz_mpoly_assert_canonical(r, ctx);

          result = fmpz_mpoly_equal(f, k, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i, j);
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

       len = n_randint(state, 15);
       len1 = n_randint(state, 15);
       len2 = n_randint(state, 15) + 1;

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_divrem_monagan_pearce(h, r, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r, ctx);
          fmpz_mpoly_remainder_test(r, g, ctx);

          fmpz_mpoly_mul_johnson(k, h, g, ctx);
	      fmpz_mpoly_add(k, k, r, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);

          result = fmpz_mpoly_equal(f, k, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check f = g*q + r for random polys\ni = %wd, j = %wd\n", i, j);
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

       len = n_randint(state, 15);
       len1 = n_randint(state, 15);
       len2 = n_randint(state, 15) + 1;

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

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

          fmpz_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          fmpz_mpoly_remainder_test(r1, g, ctx);
          fmpz_mpoly_divrem_monagan_pearce(f, r2, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r2, ctx);
          fmpz_mpoly_remainder_test(r2, g, ctx);

          result = fmpz_mpoly_equal(h, f, ctx) && fmpz_mpoly_equal(r1, r2, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of quotient with first argument\ni = %wd, j = %wd\n", i, j);
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

       len = n_randint(state, 15);
       len1 = n_randint(state, 15);
       len2 = n_randint(state, 15) + 1;

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

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

          fmpz_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          fmpz_mpoly_remainder_test(r1, g, ctx);
          fmpz_mpoly_divrem_monagan_pearce(g, r2, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r2, ctx);

          result = fmpz_mpoly_equal(h, g, ctx) && fmpz_mpoly_equal(r1, r2, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of quotient with second argument\ni = %wd, j = %wd\n", i, j);
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

       len = n_randint(state, 15);
       len1 = n_randint(state, 15);
       len2 = n_randint(state, 15) + 1;

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

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

          fmpz_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          fmpz_mpoly_remainder_test(r1, g, ctx);

          fmpz_mpoly_divrem_monagan_pearce(k, f, f, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);
          fmpz_mpoly_assert_canonical(f, ctx);
          fmpz_mpoly_remainder_test(f, g, ctx);

          result = fmpz_mpoly_equal(h, k, ctx) && fmpz_mpoly_equal(r1, f, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of remainder with first argument\ni = %wd, j = %wd\n", i, j);
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

       len = n_randint(state, 15);
       len1 = n_randint(state, 15);
       len2 = n_randint(state, 15) + 1;

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

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

          fmpz_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          fmpz_mpoly_assert_canonical(r1, ctx);
          fmpz_mpoly_remainder_test(r1, g, ctx);

          fmpz_mpoly_divrem_monagan_pearce(k, g, f, g, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);
          fmpz_mpoly_assert_canonical(g, ctx);

          result = fmpz_mpoly_equal(h, k, ctx) && fmpz_mpoly_equal(r1, g, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing of remainder with second argument\ni = %wd, j = %wd\n", i, j);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
       fmpz_mpoly_clear(k, ctx);  
       fmpz_mpoly_clear(r1, ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

