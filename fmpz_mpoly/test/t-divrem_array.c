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

    flint_printf("divrem_array....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
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

       exp_bits = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          do {
             fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(r, state, len, exp_bound, coeff_bits, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);

          ok1 = fmpz_mpoly_divrem_array(k, r, h, g, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r, g, ctx);

          result = (ok1 == 0) || (ok1 && fmpz_mpoly_equal(f, k, ctx));

          if (!result)
          {
             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k, NULL, ctx); printf("\n\n");
          
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
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

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

       exp_bits = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          do {
             fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k, state, len, exp_bound, coeff_bits, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r, f, g, ctx);
          fmpz_mpoly_test(h, ctx);
          fmpz_mpoly_test(r, ctx);
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

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");
			           
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
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

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

       exp_bits = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          do {
             fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(r1, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(r2, state, len, exp_bound, coeff_bits, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_test(h, ctx);
          fmpz_mpoly_test(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(f, r2, f, g, ctx);
          fmpz_mpoly_test(f, ctx);
          fmpz_mpoly_test(r2, ctx);
          if (ok2)
             fmpz_mpoly_remainder_test(r2, g, ctx);

          result = (ok1 == 0 || ok2 == 0) || 
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, f, ctx)
				                         && fmpz_mpoly_equal(r1, r2, ctx));

          if (!result)
          {
             printf("FAIL\n");
             printf("Aliasing test1\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r1, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r2, NULL, ctx); printf("\n\n");
          
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
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

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

       exp_bits = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          do {
             fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(r1, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(r2, state, len, exp_bound, coeff_bits, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_test(h, ctx);
          fmpz_mpoly_test(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(g, r2, f, g, ctx);
          fmpz_mpoly_test(g, ctx);
          fmpz_mpoly_test(r2, ctx);

          result = (ok1 == 0 || ok2 == 0) || 
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, g, ctx)
				                         && fmpz_mpoly_equal(r1, r2, ctx));

          if (!result)
          {
             printf("FAIL\n");
             printf("Aliasing test1\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r1, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r2, NULL, ctx); printf("\n\n");
          
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
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

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

       exp_bits = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          do {
             fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(r1, state, len, exp_bound, coeff_bits, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_test(h, ctx);
          fmpz_mpoly_test(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(k, f, f, g, ctx);
          fmpz_mpoly_test(k, ctx);
          fmpz_mpoly_test(f, ctx);
          if (ok2)
             fmpz_mpoly_remainder_test(f, g, ctx);

          result = (ok1 == 0 || ok2 == 0) || 
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, k, ctx)
				                         && fmpz_mpoly_equal(r1, f, ctx));

          if (!result)
          {
             printf("FAIL\n");
             printf("Aliasing test1\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r1, NULL, ctx); printf("\n\n");
          
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
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

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

       exp_bits = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          do {
             fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(r1, state, len, exp_bound, coeff_bits, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);

          ok1 = fmpz_mpoly_divrem_array(h, r1, f, g, ctx);
          fmpz_mpoly_test(h, ctx);
          fmpz_mpoly_test(r1, ctx);
          if (ok1)
             fmpz_mpoly_remainder_test(r1, g, ctx);

          ok2 = fmpz_mpoly_divrem_array(k, g, f, g, ctx);
          fmpz_mpoly_test(k, ctx);
          fmpz_mpoly_test(g, ctx);

          result = (ok1 == 0 || ok2 == 0) || 
                   (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, k, ctx)
				                         && fmpz_mpoly_equal(r1, g, ctx));

          if (!result)
          {
             printf("FAIL\n");
             printf("Aliasing test1\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r1, NULL, ctx); printf("\n\n");
          
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

