/*
    Copyright 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   ulong c;
   ulong B1; 
   mp_limb_t n, p;

   FLINT_TEST_INIT(state);
   

   while(1)
   {
      flint_printf("Enter number to be factored: "); fflush(stdout);
      if (!flint_scanf("%wu", &n))
      {
         flint_printf("Read failed\n");
         abort();
      }
   
      flint_printf("Enter B1: "); fflush(stdout);
      if (!flint_scanf("%wu", &B1))
      {
         flint_printf("Read failed\n");
         abort();
      }
    
      do
      {
         c = n_randint(state, n);
      } while (c <= UWORD(2));

      p = n_factor_pp1(n, B1, c);
      if (p >= 2)
         flint_printf("Factor: %wu\n", p);
      else
         flint_printf("Factor not found!\n");
   } while(1);
   
   flint_randclear(state);

   return 0;
}
