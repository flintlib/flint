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

    flint_printf("divides_monagan_pearce....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
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

       len = n_randint(state, 2);
       len1 = n_randint(state, 2);
       len2 = n_randint(state, 2) + 1;

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

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);

          ok1 = fmpz_mpoly_divides_monagan_pearce(k, h, g, ctx);
          fmpz_mpoly_test(k, ctx);

          result = (ok1 && fmpz_mpoly_equal(f, k, ctx));

          if (!result)
          {
             printf("FAIL1\n");
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
    }

    /* Check random polys don't divide */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
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

       len = n_randint(state, 20);
       len1 = n_randint(state, 20);
       len2 = n_randint(state, 20) + 1;

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

          ok1 = fmpz_mpoly_divides_monagan_pearce(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);

          if (ok1)
          {
             fmpz_mpoly_mul_johnson(k, h, g, ctx);
             fmpz_mpoly_test(k, ctx);
          }

          result = (ok1 == 0 || fmpz_mpoly_equal(f, k, ctx));

          if (!result)
          {
             printf("FAIL2\n");

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
    }

    /* Check aliasing first argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
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

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100) + 1;

       exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits1 = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits2 = n_randint(state, FLINT_BITS - 1 - 
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

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);

          ok1 = fmpz_mpoly_divides_monagan_pearce(k, h, g, ctx);
          fmpz_mpoly_test(k, ctx);
          ok2 = fmpz_mpoly_divides_monagan_pearce(h, h, g, ctx);
          fmpz_mpoly_test(h, ctx);

          result = (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, k, ctx));

          if (!result)
          {
             printf("FAIL3\n");
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
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
       fmpz_mpoly_clear(k, ctx);  
    }

    /* Check aliasing, first argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
 
       len = n_randint(state, 20);
       len1 = n_randint(state, 20);
       len2 = n_randint(state, 20) + 1;

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

          ok1 = fmpz_mpoly_divides_monagan_pearce(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);
          ok2 = fmpz_mpoly_divides_monagan_pearce(f, f, g, ctx);
          fmpz_mpoly_test(f, ctx);

          result = ((ok1 == ok2) &&  (ok1 == 0 || fmpz_mpoly_equal(f, h, ctx)));

          if (!result)
          {
             printf("FAIL4\n");

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
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
    }

    /* Check aliasing second argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
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

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100) + 1;

       exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits1 = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
       exp_bits2 = n_randint(state, FLINT_BITS - 1 - 
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

          fmpz_mpoly_mul_johnson(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);

          ok1 = fmpz_mpoly_divides_monagan_pearce(k, h, g, ctx);
          fmpz_mpoly_test(k, ctx);
          ok2 = fmpz_mpoly_divides_monagan_pearce(g, h, g, ctx);
          fmpz_mpoly_test(g, ctx);

          result = (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(g, k, ctx));

          if (!result)
          {
             printf("FAIL5\n");
             printf("Aliasing test2\n");

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
    }

    /* Check aliasing, second argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
 
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

          ok1 = fmpz_mpoly_divides_monagan_pearce(h, f, g, ctx);
          fmpz_mpoly_test(h, ctx);
          ok2 = fmpz_mpoly_divides_monagan_pearce(g, f, g, ctx);
          fmpz_mpoly_test(g, ctx);

          result = ((ok1 == ok2) &&  (ok1 == 0 || fmpz_mpoly_equal(g, h, ctx)));

          if (!result)
          {
             printf("FAIL6\n");

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

