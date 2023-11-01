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

TEST_FUNCTION_START(n_mulmod_precomp, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, b, d, r1, r2, p1, p2, dinv;
      double dpre;

      mp_limb_t bits = n_randint(state, FLINT_D_BITS) + 1;
      d = n_randtest_bits(state, bits);
      a = n_randtest(state) % d;
      b = n_randtest(state) % d;

      dpre = n_precompute_inverse(d);

      r1 = n_mulmod_precomp(a, b, d, dpre);

      umul_ppmm(p1, p2, a, b);
      dinv = n_preinvert_limb(d);
      r2 = n_ll_mod_preinv(p1, p2, d, dinv);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, b = %wu, d = %wu, dinv = %f\n", a, b, d, dpre);
         flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
