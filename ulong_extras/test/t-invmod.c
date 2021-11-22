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
   

   flint_printf("invmod....");
   fflush(stdout);
   
   for (i = 0; i < 10000 * flint_test_multiplier(); i++) 
   {
      mp_limb_t a, b, t, r, binv, ph, pl;
      
      do
      {
         a = n_randtest(state);
         b = n_randtest(state);
      } while ((a >= b) || (n_gcd(b, a) != UWORD(1)));

      t = n_invmod(a, b);
      
      binv = n_preinvert_limb(b);
      umul_ppmm(ph, pl, t, a);
      r = n_ll_mod_preinv(ph, pl, b, binv);

      result = (((r == UWORD(0)) && (b == UWORD(1))) || (r == UWORD(1)));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, b = %wu, r = %wd\n", a, b, r); 
         fflush(stdout);
         flint_abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
