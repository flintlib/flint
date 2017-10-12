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
    int i, j, w, result;
    FLINT_TEST_INIT(state);

    flint_printf("divrem_ideal_monagan_pearce....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;
       fmpz_mpoly_struct * qarr[1], * darr[1];
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

       exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits1 = n_randint(state, FLINT_BITS - 2 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits2 = n_randint(state, FLINT_BITS - 2 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

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

          qarr[0] = k;
          darr[0] = g;

          fmpz_mpoly_divrem_ideal_monagan_pearce(qarr, r, h, darr, 1, ctx);

          result = fmpz_mpoly_equal(f, k, ctx);

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

    /* Check f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, r, k1, k2;
       fmpz_mpoly_struct * g, * q;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;
       fmpz_mpoly_struct * qarr[5], * darr[5];

       num = n_randint(state, 5) + 1;

       g = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));
       q = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));

       ord = mpoly_ordering_randtest(state);
	   
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       for (w = 0; w < num; w++)
       {
          fmpz_mpoly_init(g + w, ctx);
          darr[w] = g + w;

          fmpz_mpoly_init(q + w, ctx);
          qarr[w] = q + w;
       }  

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(k1, ctx);
       fmpz_mpoly_init(k2, ctx);
       fmpz_mpoly_init(r, ctx);

       len = n_randint(state, 10);
       len1 = n_randint(state, 10);
       len2 = n_randint(state, 10) + 1;

       exp_bits = n_randint(state, 14/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 14/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 14/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          for (w = 0; w < num; w++)
          {
             do {
                fmpz_mpoly_randtest(darr[w], state, len2, exp_bound2, coeff_bits + 1, ctx);
             } while (darr[w]->length == 0);
             fmpz_mpoly_randtest(qarr[w], state, len, exp_bound, coeff_bits, ctx);
          }
          fmpz_mpoly_randtest(k1, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k2, state, len, exp_bound, coeff_bits, ctx);

          fmpz_mpoly_divrem_ideal_monagan_pearce(qarr, r, f, darr, num, ctx);

          fmpz_mpoly_zero(k2, ctx);
          for (w = 0; w < num; w++)
          {
             fmpz_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
             fmpz_mpoly_add(k2, k2, k1, ctx);
	  }
          fmpz_mpoly_add(k2, k2, r, ctx);

          result = fmpz_mpoly_equal(f, k2, ctx);

          if (!result)
          {
             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                              "coeff_bits = %ld, nvars = %ld, num = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                          len2, exp_bits2, exp_bound2, coeff_bits, nvars, num);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k2, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");

             for (w = 0; w < num; w++)
             {
                fmpz_mpoly_print_pretty(darr[w], NULL, ctx); printf("\n\n");
             }
             for (w = 0; w < num; w++)
             {
                fmpz_mpoly_print_pretty(qarr[w], NULL, ctx); printf("\n\n");
             }
			           
             flint_abort();
          }
       }

       for (w = 0; w < num; w++)
          fmpz_mpoly_clear(qarr[w], ctx);
       for (w = 0; w < num; w++)
          fmpz_mpoly_clear(darr[w], ctx);
       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(k1, ctx);  
       fmpz_mpoly_clear(k2, ctx);  
       fmpz_mpoly_clear(r, ctx);

       flint_free(g);
       flint_free(q);  
    }

    /* Check aiasing */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, r, k1, k2;
       fmpz_mpoly_struct * g, * q;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;
       fmpz_mpoly_struct * qarr[5], * darr[5];

       num = n_randint(state, 5) + 1;

       g = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));
       q = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));

       ord = mpoly_ordering_randtest(state);
	   
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       for (w = 0; w < num; w++)
       {
          fmpz_mpoly_init(g + w, ctx);
          darr[w] = g + w;

          fmpz_mpoly_init(q + w, ctx);
          qarr[w] = q + w;
       }  

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(k1, ctx);
       fmpz_mpoly_init(k2, ctx);
       fmpz_mpoly_init(r, ctx);

       len = n_randint(state, 10);
       len1 = n_randint(state, 10);
       len2 = n_randint(state, 10) + 1;

       exp_bits = n_randint(state, 14/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 14/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 14/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          for (w = 0; w < num; w++)
          {
             do {
                fmpz_mpoly_randtest(darr[w], state, len2, exp_bound2, coeff_bits + 1, ctx);
             } while (darr[w]->length == 0);
             fmpz_mpoly_randtest(qarr[w], state, len, exp_bound, coeff_bits, ctx);
          }
          fmpz_mpoly_randtest(k1, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k2, state, len, exp_bound, coeff_bits, ctx);

          fmpz_mpoly_set(r, f, ctx);

          fmpz_mpoly_divrem_ideal_monagan_pearce(qarr, f, f, darr, num, ctx);

          fmpz_mpoly_zero(k2, ctx);
          for (w = 0; w < num; w++)
          {
             fmpz_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
             fmpz_mpoly_add(k2, k2, k1, ctx);
	  }
          fmpz_mpoly_add(k2, k2, f, ctx);

          result = fmpz_mpoly_equal(r, k2, ctx);

          if (!result)
          {
             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                              "coeff_bits = %ld, nvars = %ld, num = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                          len2, exp_bits2, exp_bound2, coeff_bits, nvars, num);

             fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k2, NULL, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(r, NULL, ctx); printf("\n\n");

             for (w = 0; w < num; w++)
             {
                fmpz_mpoly_print_pretty(darr[w], NULL, ctx); printf("\n\n");
             }
             for (w = 0; w < num; w++)
             {
                fmpz_mpoly_print_pretty(qarr[w], NULL, ctx); printf("\n\n");
             }
			           
             flint_abort();
          }
       }

       for (w = 0; w < num; w++)
          fmpz_mpoly_clear(qarr[w], ctx);
       for (w = 0; w < num; w++)
          fmpz_mpoly_clear(darr[w], ctx);
       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(k1, ctx);  
       fmpz_mpoly_clear(k2, ctx);  
       fmpz_mpoly_clear(r, ctx);

       flint_free(g);
       flint_free(q);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

