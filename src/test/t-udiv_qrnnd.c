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
   

   flint_printf("udiv_qrnnd....");
   fflush(stdout);

   for (i = 0; i < 1000000; i++)
   {
      mp_limb_t d, nh, nl, q, r, ph, pl;

      do 
      {
         d = n_randtest_not_zero(state);
         nh = n_randtest(state);
      } while (nh >= d);
      nl = n_randtest(state);

      udiv_qrnnd(q, r, nh, nl, d);
      umul_ppmm(ph, pl, d, q);
      add_ssaaaa(ph, pl, ph, pl, UWORD(0), r);

      result = ((ph == nh) && (pl == nl));

      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("nh = %wu, nl = %wu, d = %wu\n", nh, nl, d); 
         flint_printf("ph = %wu, pl = %wu\n", ph, pl);
         fflush(stdout);
         flint_abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
