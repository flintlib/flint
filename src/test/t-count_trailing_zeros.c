/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   FLINT_TEST_INIT(state);
   

   flint_printf("count_trailing_zeros....");
   fflush(stdout);

   for (i = 0; i < 1000000; i++)
   {
      mp_limb_t n;
      unsigned int count = 0;

      n = n_randtest(state);

      if (n != 0)
         count_trailing_zeros(count, n);

      result = ((n == UWORD(0)) || (((n >> count) & UWORD(1)) && (l_shift(n, FLINT_BITS-count) == UWORD(0))));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wu, count = %u\n", n, count); 
         fflush(stdout);
         flint_abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
