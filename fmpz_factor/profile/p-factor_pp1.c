/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright 2012 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_factor.h"

int main(void)
{
   fmpz_t n, p;
   ulong c;
   ulong B1;
   
   fmpz_init(n);
   fmpz_init(p);

   flint_rand_t state;
   flint_randinit(state);

   while(1)
   {
      printf("Enter number to be factored: "); fflush(stdout);
      if (!fmpz_read(n))
      {
         printf("Read failed\n");
         abort();
      }
   
      printf("Enter B1: "); fflush(stdout);
      if (!scanf("%lu", &B1))
      {
         printf("Read failed\n");
         abort();
      }
    
      do
      {
         c = n_randlimb(state);
      } while (c <= 2UL);

      if (fmpz_factor_pp1(p, n, B1, c))
      {
         printf("Factor: ");
         fmpz_print(p);
         printf("\n");
      } else
         printf("Factor not found!\n");
   } while(1);
   
   flint_randclear(state);

   fmpz_clear(n);
   fmpz_clear(p);

   return 0;
}