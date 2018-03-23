/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result, ok1, ok2;
    FLINT_TEST_INIT(state);

    flint_printf("mul_array....");
    fflush(stdout);

    /* Check f*g = g*f */
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
       len2 = n_randint(state, 100);

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_mul_array(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          ok2 = fmpz_mpoly_mul_array(k, g, f, ctx);
          fmpz_mpoly_assert_canonical(k, ctx);

          result = (ok1 == 0 && ok2 == 0) || fmpz_mpoly_equal(h, k, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check f*g = g*f\ni = %wd, j = %wd\n", i, j);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
       fmpz_mpoly_clear(k, ctx);  
    }

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k1, k2, t1, t2;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(t1, ctx);
       fmpz_mpoly_init(t2, ctx);
       fmpz_mpoly_init(k1, ctx);
       fmpz_mpoly_init(k2, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(k1, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k2, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
          fmpz_mpoly_randtest_bound(h, state, len2, coeff_bits, exp_bound2, ctx);

          fmpz_mpoly_add(t1, g, h, ctx);
          ok1 = fmpz_mpoly_mul_array(k1, f, t1, ctx);
          fmpz_mpoly_assert_canonical(k1, ctx);

          ok2 = fmpz_mpoly_mul_array(t1, f, g, ctx);
          fmpz_mpoly_assert_canonical(t1, ctx);
          if (ok2)
          {
             ok2 = fmpz_mpoly_mul_array(t2, f, h, ctx);
             fmpz_mpoly_assert_canonical(t2, ctx);
          }
          if (ok2)
          {
             fmpz_mpoly_add(k2, t1, t2, ctx);
             fmpz_mpoly_assert_canonical(k2, ctx);
          }

          result = (ok1 == 0 || ok2 == 0) || fmpz_mpoly_equal(k1, k2, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i, j);
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
    }

    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_mul_array(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok2 = fmpz_mpoly_mul_array(f, f, g, ctx);
          fmpz_mpoly_assert_canonical(f, ctx);

          result = (ok1 == 0 && ok2 == 0) || fmpz_mpoly_equal(h, f, ctx);

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
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       exp_bound =  n_randint(state, 800/nvars/nvars) + 1;
       exp_bound1 = n_randint(state, 800/nvars/nvars) + 1;
       exp_bound2 = n_randint(state, 800/nvars/nvars) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_mul_array(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);

          ok2 = fmpz_mpoly_mul_array(g, f, g, ctx);
          fmpz_mpoly_assert_canonical(g, ctx);

          result = (ok1 == 0 && ok2 == 0) || fmpz_mpoly_equal(h, g, ctx);

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
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

