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

TEST_FUNCTION_START(n_pow, state)
{
   int i, result;

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test a^e1 * a^e2 = a^(e1 + e2) */
   {
      mp_limb_t exp1, exp2, n, bits, r1, r2;

      bits = n_randint(state, 55) + 10;
      exp1 = n_randint(state, 5);
      exp2 = n_randint(state, 5);

      if ((exp1 == WORD(0)) && (exp2 == WORD(0))) bits = n_randint(state, 64) + 1;
      else bits /= (exp1 + exp2);

      n = n_randtest_bits(state, bits);

      r1 = n_pow(n, exp1)*n_pow(n, exp2);
      r2 = n_pow(n, exp1 + exp2);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wu, exp1 = %wu, exp2 = %wu, r1 = %wu, r2 = %wu\n", n, exp1, exp2, r1, r2);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
