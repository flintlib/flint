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
   
   flint_printf("gcd_full....");
   fflush(stdout);
   
   

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) 
   {
      mp_limb_t a, b, c, bits1, bits2, bits3;
      
      bits1 = n_randint(state, FLINT_BITS-1) + 1;
      bits2 = n_randint(state, bits1) + 1;
      bits3 = n_randint(state, FLINT_BITS - bits1) + 1;

      do
      {
         a = n_randtest_bits(state, bits1);
         b = n_randtest_bits(state, bits2);
      } while ((n_gcd_full(a, b) != UWORD(1)));

      c = n_randtest_bits(state, bits3);

      result = (n_gcd_full(a*c, b*c) == c);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, b = %wu, c = %wu\n", a, b, c); 
         abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
