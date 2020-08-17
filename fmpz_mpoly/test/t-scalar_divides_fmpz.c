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
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_divides_fmpz....");
    fflush(stdout);

    /* Check (f*a)/a = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       fmpz_t c;
       slong nvars, len, coeff_bits, exp_bits;

       fmpz_init(c);

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 10; j++)
       {
          fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
          fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
          fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

          fmpz_randtest_not_zero(c, state, n_randint(state, 200) + 1);

          fmpz_mpoly_scalar_mul_fmpz(g, f, c, ctx);
          result = fmpz_mpoly_scalar_divides_fmpz(h, g, c, ctx);
          result = result && fmpz_mpoly_equal(h, f, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check (f*a)/a = f\ni = %wd, j = %wd\n", i,j);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  

       fmpz_clear(c);      
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       fmpz_t c;
       slong nvars, len, coeff_bits, exp_bits;

       fmpz_init(c);

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 10; j++)
       {
          fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
          fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
       
          fmpz_randtest(c, state, n_randint(state, 200));

          fmpz_mpoly_scalar_mul_fmpz(f, f, c, ctx);
          
          fmpz_mpoly_set(g, f, ctx);

          result = fmpz_mpoly_scalar_divides_fmpz(h, f, c, ctx);
          result = result && fmpz_mpoly_scalar_divides_fmpz(g, g, c, ctx);
          result = result && fmpz_mpoly_equal(g, h, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing\ni = %wd, j = %wd\n", i,j);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  

       fmpz_clear(c);      
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

