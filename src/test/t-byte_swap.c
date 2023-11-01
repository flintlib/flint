/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

ulong byte_swap_naive(ulong n)
{
   ulong r = 0;
   slong i;

   for (i = 0; i < sizeof(ulong); i++)
   {
      r <<= 8;
      r |= (n & 0xFF);
      n >>= 8;
   }

   return r;
}


TEST_FUNCTION_START(byte_swap, state)
{
   int i, result;

   /* byte_swap(byte_swap(n)) == n */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      ulong n, r1, r2;

      n = n_randtest(state);

      r1 = n;

      r2 = n;
      byte_swap(r2);
      byte_swap(r2);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("byte_swap(byte_swap(n)) != n\n");
         flint_printf("n = %wx, r1 = %wx, r2 = %wx\n", n, r1, r2);
         fflush(stdout);
         flint_abort();
      }
   }

   /* byte_swap(n) == byte_swap_naive(n) */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      ulong n, r1, r2;

      n = n_randtest(state);

      r1 = n;
      byte_swap(r1);

      r2 = byte_swap_naive(n);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("byte_swap(n) != byte_swap_naive(n)\n");
         flint_printf("n = %wx, r1 = %wx, r2 = %wx\n", n, r1, r2);
         fflush(stdout);
         flint_abort();
      }
   }

   TEST_FUNCTION_END(state);
}
