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

void fmpz_mpoly_pow_naive(fmpz_mpoly_t res, fmpz_mpoly_t f,
                                                 slong n, fmpz_mpoly_ctx_t ctx)
{
   if (f->length == 0)
      fmpz_mpoly_zero(res, ctx);
   else if (n == 0)
      fmpz_mpoly_set_ui(res, 1, ctx);
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

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("pow_fps....");
    fflush(stdout);

    /* Check pow_fps against pow_naive */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong nvars, len, len1, exp_bound, exp_bound1;
       slong coeff_bits, exp_bits, exp_bits1, s;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       s = n_randint(state, 12);

       len = n_randint(state, 10);
       len1 = n_randint(state, 10);

       exp_bits = n_randint(state, FLINT_BITS - 1 -
                         mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits1 = n_randint(state, FLINT_BITS - 1 - FLINT_BIT_COUNT(s) -
                         mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          fmpz_mpoly_randtest(g, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
             
          fmpz_mpoly_pow_fps(g, f, s, ctx);
          fmpz_mpoly_test(g, ctx);

          fmpz_mpoly_pow_naive(h, f, s, ctx);
          fmpz_mpoly_test(h, ctx);

          result = fmpz_mpoly_equal(g, h, ctx);

          if (!result)
          {
             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                                                            coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g;
       ordering_t ord;
       slong nvars, len, len1, exp_bound, exp_bound1;
       slong coeff_bits, exp_bits, exp_bits1, s;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);

       s = n_randint(state, 12);

       len = n_randint(state, 10);
       len1 = n_randint(state, 10);

       exp_bits = n_randint(state, FLINT_BITS - 1 -
                         mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits1 = n_randint(state, FLINT_BITS - 1 - FLINT_BIT_COUNT(s) -
                         mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          fmpz_mpoly_randtest(g, state, len, exp_bound, coeff_bits, ctx);

          fmpz_mpoly_pow_fps(g, f, s, ctx);
          fmpz_mpoly_test(g, ctx);

          fmpz_mpoly_pow_fps(f, f, s, ctx);
          fmpz_mpoly_test(f, ctx);

          result = fmpz_mpoly_equal(f, g, ctx);

          if (!result)
          {
             printf("FAIL\n");
             printf("Aliasing test1\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                                                            coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

