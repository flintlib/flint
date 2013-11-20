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
   FLINT_TEST_INIT(state);
   
   flint_printf("factor_power235....");
   fflush(stdout);
 
   

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares */
   {
      mp_limb_t factor, exp, n1, n2, bits;
      
      bits = n_randint(state, FLINT_BITS/2) + 1;
      n1 = n_randtest_bits(state, bits);
      factor = n_factor_power235(&exp, n1*n1);

      n2 = n_pow(factor, exp);

      result = (n1*n1 == n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("factor = %wu, exp = %wu\n", factor, exp); 
         abort();
      }
   }
   
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random cubes */
   {
      mp_limb_t factor, exp, n1, n2, bits;
      
      bits = n_randint(state, FLINT_BITS/3) + 1;
      n1 = n_randtest_bits(state, bits);
      factor = n_factor_power235(&exp, n1*n1*n1);

      n2 = n_pow(factor, exp);

      result = (n1*n1*n1 == n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("factor = %wu, exp = %wu\n", factor, exp); 
         abort();
      }
   }
   
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random fifth powers */
   {
      mp_limb_t factor, exp, n1, n2, bits;
      
      bits = n_randint(state, FLINT_BITS/5) + 1;
      n1 = n_randtest_bits(state, bits);
      factor = n_factor_power235(&exp, n1*n1*n1*n1*n1);

      n2 = n_pow(factor, exp);

      result = (n1*n1*n1*n1*n1 == n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("factor = %wu, exp = %wu\n", factor, exp); 
         abort();
      }
   }
   
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test non 235-powers */
   {
      mp_limb_t exp, n1;
      
      do
      {
         n1 = n_randtest(state);
      } while (n_is_perfect_power235(n1));
      
      result = (!n_factor_power235(&exp, n1));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n1 = %wu, exp = %wu\n", n1, exp); 
         abort();
      }
   }
   
   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
