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

TEST_FUNCTION_START(n_is_oddprime_binary, state)
{
   int i, result;
   slong cutoff = 100000;

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that primes pass the test */
   {
      mp_limb_t d;
      mpz_t d_m;

      mpz_init(d_m);

      do
      {
         d = n_randint(state, cutoff) | 1;
         if (d == UWORD(1)) d += 16; /* algorithm requires d >= 17 */
         flint_mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = flint_mpz_get_ui(d_m);
      } while (d > cutoff);

      result = n_is_oddprime_binary(d);

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
         d = (n_randint(state, cutoff) + 16) | 1;
         flint_mpz_set_ui(d_m, d);
      } while ((mpz_probab_prime_p(d_m, 12)) || (d > cutoff));

      result = !n_is_oddprime_binary(d);

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
