/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_factor_power235, state)
{
   int i, result;

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares */
   {
      mp_limb_t factor, exp, n1, n2, bits;

      bits = n_randint(state, FLINT_BITS/2) + 1;
      n1 = n_randtest_bits(state, bits);
      factor = n_factor_power235(&exp, n1*n1);

      n2 = n_pow(factor, exp);

      result = (n1*n1 == n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("factor = %wu, exp = %wu\n", factor, exp);
         fflush(stdout);
         flint_abort();
      }
   }

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random cubes */
   {
      mp_limb_t factor, exp, n1, n2, bits;

      bits = n_randint(state, FLINT_BITS/3) + 1;
      n1 = n_randtest_bits(state, bits);
      factor = n_factor_power235(&exp, n1*n1*n1);

      n2 = n_pow(factor, exp);

      result = (n1*n1*n1 == n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("factor = %wu, exp = %wu\n", factor, exp);
         fflush(stdout);
         flint_abort();
      }
   }

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random fifth powers */
   {
      mp_limb_t factor, exp, n1, n2, bits;

      bits = n_randint(state, FLINT_BITS/5) + 1;
      n1 = n_randtest_bits(state, bits);
      factor = n_factor_power235(&exp, n1*n1*n1*n1*n1);

      n2 = n_pow(factor, exp);

      result = (n1*n1*n1*n1*n1 == n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("factor = %wu, exp = %wu\n", factor, exp);
         fflush(stdout);
         flint_abort();
      }
   }

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test non 235-powers */
   {
      mp_limb_t exp, n1;

      do
      {
         n1 = n_randtest(state);
      } while (n_is_perfect_power235(n1));

      result = (!n_factor_power235(&exp, n1));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n1 = %wu, exp = %wu\n", n1, exp);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
