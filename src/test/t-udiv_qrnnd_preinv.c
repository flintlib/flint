/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(udiv_qrnnd_preinv, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t d, dinv, nh, nl, q1, r1, q2, r2, norm;

      do
      {
         d = n_randtest_not_zero(state);
         nh = n_randtest(state);
         norm = flint_clz(d);
         d <<= norm;
      } while (nh >= d);
      nl = n_randtest(state);

      dinv = n_preinvert_limb_prenorm(d);

      udiv_qrnnd_preinv(q1, r1, nh, nl, d, dinv);
      udiv_qrnnd(q2, r2, nh, nl, d);

      result = ((q1 == q2) && (r1 == r2));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("nh = %wu, nl = %wu, d = %wu, dinv = %wu\n", nh, nl, d, dinv);
         flint_printf("q1 = %wu, q2 = %wu, r1 = %wu, r2 = %wu\n", q1, q2, r1, r2);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
