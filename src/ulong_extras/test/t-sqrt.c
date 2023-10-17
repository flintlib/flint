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

TEST_FUNCTION_START(n_sqrt, state)
{
   int i, result;

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, s1, s2;
      mpz_t a_m, s2_m;

      mpz_init(a_m);
      mpz_init(s2_m);

      a = n_randtest(state);

      s1 = n_sqrt(a);

      flint_mpz_set_ui(a_m, a);
      mpz_sqrt(s2_m, a_m);
      s2 = flint_mpz_get_ui(s2_m);

      result = (s1 == s2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, s1 = %wd, s2 = %wu\n", a, s1, s2);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(a_m);
      mpz_clear(s2_m);
   }

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, s1, s2, bits;
      mpz_t a_m, s2_m;

      mpz_init(a_m);
      mpz_init(s2_m);

      bits = n_randint(state, FLINT_BITS/2 + 1);
      a = n_randtest_bits(state, bits);
      a = a*a;
      a += (n_randint(state, 100) - 50);
      s1 = n_sqrt(a);

      flint_mpz_set_ui(a_m, a);
      mpz_sqrt(s2_m, a_m);
      s2 = flint_mpz_get_ui(s2_m);

      result = (s1 == s2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, s1 = %wd, s2 = %wu\n", a, s1, s2);
         fflush(stdout);
         flint_abort();
      }

      mpz_clear(a_m);
      mpz_clear(s2_m);
   }

   TEST_FUNCTION_END(state);
}
