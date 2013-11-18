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
   int i, j, result;
   FLINT_TEST_INIT(state);
   

   flint_printf("factor_trial....");
   fflush(stdout);
 
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      mp_limb_t n1, n2;
      n_factor_t factors;

      n_factor_init(&factors);

      n1 = n_randtest_not_zero(state);
      n2 = n_factor_trial(&factors, n1, UWORD(10000));
      
      for (j = 0; j < factors.num; j++)
      {
         n2 *= n_pow(factors.p[j], factors.exp[j]);
      }

      result = (n1 == n2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n1 = %wu, n2 = %wu\n", n1, n2); 
         abort();
      }
   }
   
   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
