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
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_mul_fmpz....");
    fflush(stdout);

    /* Check (f*a)*b = f*(a*b) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
       ordering_t ord;
       fmpz_t a, b, c;
       slong nvars, len, exp_bound, coeff_bits, exp_bits;

       fmpz_init(a);
       fmpz_init(b);
       fmpz_init(c);

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, FLINT_BITS -
                     mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars) - 1) + 1;
       exp_bound = n_randbits(state, exp_bits);
       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 10; j++)
       {
          fmpz_mpoly_randtest(f, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(g, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k, state, len, exp_bound, coeff_bits, ctx);

          fmpz_randtest(a, state, n_randint(state, 200));
          fmpz_randtest(b, state, n_randint(state, 200));
          fmpz_mul(c, a, b);

          fmpz_mpoly_scalar_mul_fmpz(g, f, a, ctx);
          fmpz_mpoly_scalar_mul_fmpz(h, g, b, ctx);

          fmpz_mpoly_scalar_mul_fmpz(k, f, c, ctx);

          result = fmpz_mpoly_equal(h, k, ctx);

          if (!result)
          {
             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                                  len, exp_bits, exp_bound, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k, NULL, ctx); printf("\n\n");
          
             printf("a = "); fmpz_print(a); printf("\n");
             printf("b = "); fmpz_print(b); printf("\n");
             printf("c = "); fmpz_print(c); printf("\n");

             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
       fmpz_mpoly_clear(k, ctx);  

       fmpz_clear(a);      
       fmpz_clear(b);      
       fmpz_clear(c);      
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       fmpz_t c;
       slong nvars, len, exp_bound, coeff_bits, exp_bits;

       fmpz_init(c);

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, FLINT_BITS -
                     mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars) - 1) + 1;
       exp_bound = n_randbits(state, exp_bits);
       coeff_bits = n_randint(state, 200);
       
       for (j = 0; j < 10; j++)
       {
          fmpz_mpoly_randtest(f, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);

          fmpz_randtest(c, state, n_randint(state, 200));

          fmpz_mpoly_set(g, f, ctx);

          fmpz_mpoly_scalar_mul_fmpz(h, f, c, ctx);
          
          fmpz_mpoly_scalar_mul_fmpz(g, g, c, ctx);

          result = fmpz_mpoly_equal(g, h, ctx);

          if (!result)
          {
             printf("FAIL\n");
             
             printf("Aliasing test\n");
             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                                  len, exp_bits, exp_bound, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
          
             printf("c = "); fmpz_print(c); printf("\n");

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

