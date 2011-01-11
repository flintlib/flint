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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   flint_rand_t state;
   
   printf("gcd....");
   fflush(stdout);
   
   flint_randinit(state);

   for (i = 0; i < 100000; i++) 
   {
      mp_limb_t a, b, c, bits1, bits2, bits3;
      
      bits1 = n_randint(state, FLINT_BITS-1) + 1;
      bits2 = n_randint(state, bits1) + 1;
      bits3 = n_randint(state, FLINT_BITS - bits1) + 1;

      do
      {
         a = n_randbits(bits1, state);
         b = n_randbits(bits2, state);
      } while ((n_gcd(a, b) != 1UL) || (b > a));

      c = n_randbits(bits3, state);

      result = (n_gcd(a*c, b*c) == c);
      if (!result)
      {
         printf("FAIL:\n");
         printf("a = %lu, b = %lu, c = %lu\n", a, b, c); 
         abort();
      }
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
