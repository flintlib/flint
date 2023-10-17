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

TEST_FUNCTION_START(n_factor_one_line, state)
{
   int i, result;
   slong num_iter;
   ulong count = UWORD(0);

   num_iter = 500 * FLINT_MAX(1, flint_test_multiplier());

   for (i = 0; i < num_iter; i++) /* Test random numbers */
   {
      mp_limb_t n1, n2, bits;

      do
      {
#if FLINT64
         bits = n_randint(state, 44);
#else
         bits = n_randint(state, 20);
#endif
         n1 = n_randtest_bits(state, bits + 1);
      } while (n_is_prime(n1) || (n1 == UWORD(1)));

      n2 = n_factor_one_line(n1, 50000);

      if (n2)
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

   if (count < 0.9 * num_iter)
   {
      flint_printf("FAIL:\n");
      flint_printf("Only %wu numbers factored\n", count);
      fflush(stdout);
      flint_abort();
   }

   TEST_FUNCTION_END(state);
}
