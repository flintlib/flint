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
   int i, j, result;
   FLINT_TEST_INIT(state);
   

   flint_printf("umul_ppmm....");
   fflush(stdout);

   for (i = 0; i < 1000000; i++)
   {
      mp_limb_t ph1, pl1, ph2, pl2, pl2old, m1, m2, bit;

      m1 = n_randtest(state);
      m2 = n_randtest(state);
      
      umul_ppmm(ph1, pl1, m1, m2);
      
      pl2old = UWORD(0);
      pl2 = UWORD(0);
      ph2 = UWORD(0);
      bit = UWORD(1);
      for (j = 0; j < FLINT_BITS; j++)
      {
         if (m2 & bit)
         {
            pl2 += (m1 << j);
            ph2 += (pl2 < pl2old);
            ph2 += r_shift(m1, FLINT_BITS - j);
            pl2old = pl2;
         }
         bit <<= 1;
      }

      result = ((ph2 == ph1) && (pl2 == pl1));

      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("m1 = %wu, m2 = %wu\n", m1, m2); 
         flint_printf("ph2 = %wu, ph1 = %wu, pl2 = %wu, pl1 = %wu\n", ph2, ph1, pl2, pl1);
         flint_abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
