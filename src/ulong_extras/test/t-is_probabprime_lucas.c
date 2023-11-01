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

TEST_FUNCTION_START(n_is_probabprime_lucas, state)
{
   int i, result;
   ulong count = UWORD(0);
   mp_limb_t d;
   mpz_t d_m;
   slong test_multiplier;

   test_multiplier = FLINT_MAX(1, flint_test_multiplier());

   for (i = 0; i < 10000 * test_multiplier; i++) /* Test that primes pass the test */
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest_not_zero(state);
         flint_mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = flint_mpz_get_ui(d_m);
      } while (mpz_size(d_m) > 1);

      result = n_is_probabprime_lucas(d);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("d = %wu is declared composite\n", d);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(d_m);
   }

   for (i = 0; i < 10000 * test_multiplier; i++) /* Test that not too many composites pass */
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest(state);
         flint_mpz_set_ui(d_m, d);
      } while (mpz_probab_prime_p(d_m, 12));

      if (n_is_probabprime_lucas(d) == 1) count++;

      mpz_clear(d_m);
   }

   result = (count < 20 * test_multiplier);
   if (!result)
   {
      flint_printf("FAIL:\n");
      flint_printf("%wu composites declared prime\n", count);
      fflush(stdout);
      flint_abort();
   }

   TEST_FUNCTION_END(state);
}
