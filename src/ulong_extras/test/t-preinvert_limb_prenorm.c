/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

#define invert_limb_naive(ninv, n)                    \
   do {                                               \
      mp_limb_t dummy;                                \
      udiv_qrnnd (ninv, dummy, ~(n), ~(WORD(0)), n);  \
   } while (0)

TEST_FUNCTION_START(n_preinvert_limb_prenorm, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t n, ninv1, ninv2;

      n = n_randtest(state);
      n |= (UWORD(1) << (FLINT_BITS - 1));

      ninv1 = n_preinvert_limb_prenorm(n);
      invert_limb_naive(ninv2, n);

      result = (ninv1 == ninv2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wx, ninv1 = %wx, ninv2 = %wx\n", n, ninv1, ninv2);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state)
}
