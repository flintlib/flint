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
   

   flint_printf("mod2_precomp....");
   fflush(stdout);

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t d, n, r1, r2;
      double dpre;

      d = n_randtest_not_zero(state);
      
      n = n_randtest(state);
      
      dpre = n_precompute_inverse(d);

      r1 = n_mod2_precomp(n, d, dpre);
      r2 = n%d;

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wu, d = %wu, dpre = %f\n", n, d, dpre); 
         flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
         abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
