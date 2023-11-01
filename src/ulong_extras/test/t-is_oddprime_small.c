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

TEST_FUNCTION_START(n_is_oddprime_small, state)
{
   int i, result;

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that primes pass the test */
   {
      mp_limb_t d;
      mpz_t d_m;

      mpz_init(d_m);

      do
      {
         d = n_randint(state, FLINT_ODDPRIME_SMALL_CUTOFF) | 1;
         flint_mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = flint_mpz_get_ui(d_m);
      } while (d > FLINT_ODDPRIME_SMALL_CUTOFF);

      result = n_is_oddprime_small(d);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("d = %wu is declared composite\n", d);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(d_m);
   }

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that not too many composites pass */
   {
      mp_limb_t d;
      mpz_t d_m;

      mpz_init(d_m);

      do
      {
         d = n_randint(state, FLINT_ODDPRIME_SMALL_CUTOFF) | 1;
         flint_mpz_set_ui(d_m, d);
      } while ((mpz_probab_prime_p(d_m, 12)) || (d > FLINT_ODDPRIME_SMALL_CUTOFF));

      result = !n_is_oddprime_small(d);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("d = %wu is declared prime\n", d);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(d_m);
   }

   TEST_FUNCTION_END(state);
}
