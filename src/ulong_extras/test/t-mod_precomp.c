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

TEST_FUNCTION_START(n_mod_precomp, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t bits, d, n, r1, r2;
      double dpre;

      bits = n_randint(state, FLINT_D_BITS) + 1;
      d = n_randtest_bits(state, bits);
      if (bits <= (FLINT_BITS/2)) n = n_randtest(state) % (d*d);
      else n = n_randtest(state);

      /* must have n < 2^(FLINT_BITS - 1) */
      if (FLINT_BIT_COUNT(n) == FLINT_BITS)
         n >>= 1;

      dpre = n_precompute_inverse(d);

      r1 = n_mod_precomp(n, d, dpre);
      r2 = n%d;

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wu, d = %wu, dinv = %g\n", n, d, dpre);
         flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
