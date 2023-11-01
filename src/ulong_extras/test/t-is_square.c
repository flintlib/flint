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

TEST_FUNCTION_START(n_is_square, state)
{
   int i, result;

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that non-squares pass */
   {
      mp_limb_t a, s, bits;

      bits = n_randint(state, FLINT_BITS/2) + 1;
      a = n_randtest_bits(state, bits);
      s = a*a + n_randtest(state) % (2*a) + 1;

      result = !n_is_square(s);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("s = %wu is declared square\n", s);
         fflush(stdout);
         flint_abort();
      }
   }

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that squares pass */
   {
      mp_limb_t a, s, bits;

      bits = n_randint(state, FLINT_BITS/2);
      a = n_randtest_bits(state, bits);
      s = a*a;

      result = n_is_square(s);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("s = %wu is declared square\n", s);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
