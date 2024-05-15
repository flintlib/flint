/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

#define l_shift(in, shift) \
    ((shift == FLINT_BITS) ? WORD(0) : ((in) << (shift)))

TEST_FUNCTION_START(flint_ctz, state)
{
   int i, result;

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      ulong n;
      unsigned int count;

      n = n_randtest(state);

      if (n == 0)
         continue;

      count = flint_ctz(n);

      result = ((n >> count) & UWORD(1)) && (l_shift(n, FLINT_BITS-count) == UWORD(0));
      if (!result)
            TEST_FUNCTION_FAIL("n = %wu, count = %u\n", n, count);
   }

   TEST_FUNCTION_END(state);
}
