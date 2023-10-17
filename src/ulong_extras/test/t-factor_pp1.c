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

TEST_FUNCTION_START(n_factor_pp1, state)
{
   int i, j, result;
   ulong count = UWORD(0);

   for (i = 0; i < 300 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      mp_limb_t n1, n2;

      do
      {
         n1 = n_randtest_bits(state, n_randint(state, FLINT_BITS) + 1);
      } while (n_is_prime(n1) || (n1 < UWORD(2)));

      for (j = 0; j < 20; j++)
      {
         n2 = n_factor_pp1(n1, 1000000, n_randint(state, n1 - 3) + 3);
         if (n2 > 1) break;
      }

      if (n2 > 1)
      {
         count++;
         result = ((n1%n2) == UWORD(0));
         if (!result)
         {
            flint_printf("FAIL:\n");
            flint_printf("n1 = %wu, n2 = %wu\n", n1, n2);
            fflush(stdout);
            flint_abort();
         }
      }
   }

   if (count < 295 * flint_test_multiplier())
   {
      flint_printf("FAIL:\n");
      flint_printf("Only %wu numbers factored\n", count);
      fflush(stdout);
      flint_abort();
   }

   TEST_FUNCTION_END(state);
}
