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
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   flint_rand_t state;
   
   printf("pow....");
   fflush(stdout);
 
   flint_randinit(state);

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test a^e1 * a^e2 = a^(e1 + e2) */
   {
      mp_limb_t exp1, exp2, n, bits, r1, r2;
      
      bits = n_randint(state, 55) + 10;
      exp1 = n_randint(state, 5);
      exp2 = n_randint(state, 5);
      
      if ((exp1 == 0L) && (exp2 == 0L)) bits = n_randint(state, 64) + 1;
      else bits /= (exp1 + exp2);
      
      n = n_randtest_bits(state, bits);
      
      r1 = n_pow(n, exp1)*n_pow(n, exp2);
      r2 = n_pow(n, exp1 + exp2);

      result = (r1 == r2);
      if (!result)
      {
         printf("FAIL:\n");
         printf("n = %lu, exp1 = %lu, exp2 = %lu, r1 = %lu, r2 = %lu\n", n, exp1, exp2, r1, r2); 
         abort();
      }
   }
   
   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
