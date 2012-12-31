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
   long iters;
   
   fmpz_init(n);
   fmpz_init(p);

   printf("Enter number to be factored: "); fflush(stdout);
   if (!fmpz_read(n))
   {
      printf("Read failed\n");
      abort();
   }
   
   printf("Enter a number of iterations: "); fflush(stdout);
   if (!scanf("%ld", &iters))
   {
      printf("Read failed\n");
      abort();
   }
    
   /* find prime such that n is a square mod p (or p divides n) */
   if (fmpz_is_even(n))
   {
      printf("Factor: 2\n");
      goto cleanup;
   }
   
   if (fmpz_factor_pp1(p, n, iters))
   {
      printf("Factor: ");
      fmpz_print(p);
      printf("\n");
   } else
      printf("Factor not found!\n");
   
cleanup:

   fmpz_clear(n);
   fmpz_clear(p);

   return 0;
}