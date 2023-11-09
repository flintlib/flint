/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_remove2_precomp, state)
{
   int i, result;
   const mp_limb_t * primes;
   const double * inverses;

   primes = n_primes_arr_readonly(10000);
   inverses = n_prime_inverses_arr_readonly(10000);

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      mp_limb_t n1, n2, orig_n;
      mpz_t d_n2, d_n1, d_p;
      int exp1, exp2;
      ulong j;

      mpz_init(d_n1);
      mpz_init(d_n2);
      mpz_init(d_p);

      n1 = n_randtest_not_zero(state);
      orig_n = n1;

      for (j = 0; j < FLINT_NUM_PRIMES_SMALL/10; j++)
      {
         flint_mpz_set_ui(d_n1, n1);
         flint_mpz_set_ui(d_p, flint_primes_small[j]);
         exp1 = n_remove2_precomp(&n1, primes[j], inverses[j]);
         exp2 = mpz_remove(d_n2, d_n1, d_p);
         n2 = flint_mpz_get_ui(d_n2);

         result = ((exp1 == exp2) && (n1 == n2));
         if (!result)
         {
            flint_printf("FAIL\n");
            flint_printf("n = %wu, exp1 = %d, exp2 = %d, n1 = %wu, n2 = %wu, p = %d\n", orig_n, exp1, exp2, n1, n2, flint_primes_small[j]);
            fflush(stdout);
            flint_abort();
         }
      }

      mpz_clear(d_n1);
      mpz_clear(d_n2);
      mpz_clear(d_p);
   }

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test perfect powers */
   {
      mp_limb_t n1, n2, orig_n, base;
      mpz_t d_n2, d_n1, d_p;
      int exp1, exp2, exp;
      ulong j;

      mpz_init(d_n1);
      mpz_init(d_n2);
      mpz_init(d_p);

      base = n_randtest_not_zero(state);
      n1 = base;
      exp = n_randint(state, FLINT_BITS/FLINT_BIT_COUNT(n1)) + 1;
      n1 = n_pow(base, exp);

      orig_n = n1;

      for (j = 0; j < FLINT_NUM_PRIMES_SMALL/10; j++)
      {
         flint_mpz_set_ui(d_n1, n1);
         flint_mpz_set_ui(d_p, flint_primes_small[j]);
         exp1 = n_remove2_precomp(&n1, primes[j], inverses[j]);
         exp2 = mpz_remove(d_n2, d_n1, d_p);
         n2 = flint_mpz_get_ui(d_n2);

         result = ((exp1 == exp2) && (n1 == n2));
         if (!result)
         {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu, exp1 = %d, exp2 = %d, n1 = %wu, n2 = %wu, p = %d\n", orig_n, exp1, exp2, n1, n2, flint_primes_small[j]);
            fflush(stdout);
            flint_abort();
         }
      }

      mpz_clear(d_n1);
      mpz_clear(d_n2);
      mpz_clear(d_p);
   }

   TEST_FUNCTION_END(state);
}
