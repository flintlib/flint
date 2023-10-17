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

TEST_FUNCTION_START(n_factor_partial, state)
{
   int i, j, result;

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      mp_limb_t n1, n2, prod, limit;
      n_factor_t factors;

      n_factor_init(&factors);

      n1 = n_randtest_not_zero(state);
      limit = n_sqrt(n1);
      n2 = n_factor_partial(&factors, n1, limit, 0);

      prod = 1;
      for (j = 0; j < factors.num; j++)
      {
         prod *= n_pow(factors.p[j], factors.exp[j]);
      }

      result = ((n1 == n2*prod) && ((prod > limit) || (n1 == 1)));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n1 = %wu, n2 = %wu\n", n1, n2);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
