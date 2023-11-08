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

TEST_FUNCTION_START(flint_ctz, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t n;
      unsigned int count = 0;

      n = n_randtest(state);

      if (n != 0)
         count = flint_ctz(n);

      result = ((n == UWORD(0)) || (((n >> count) & UWORD(1)) && (l_shift(n, FLINT_BITS-count) == UWORD(0))));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wu, count = %u\n", n, count);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
