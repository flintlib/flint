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
   ulong bits;
   mp_limb_t d;
   flint_rand_t state;
   
   printf("is_perfect_power235....");
   fflush(stdout);
   
   flint_randinit(state);

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that square pass the test */
   {
      bits = n_randint(state, FLINT_BITS/2) + 1;
      d = n_randtest_bits(state, bits);

      result = n_is_perfect_power235(n_pow(d, 2));
      if (!result)
      {
         printf("FAIL:\n");
         printf("d^2 = %lu is declared not a perfect power\n", d*d); 
         abort();
      }

   }
         
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that cubes pass the test */
   {
      bits = n_randint(state, FLINT_BITS/3) + 1;
      d = n_randtest_bits(state, bits);

      result = n_is_perfect_power235(n_pow(d, 3));
      
      if (!result)
      {
         printf("FAIL:\n");
         printf("d^3 = %lu is declared not a perfect power\n", d*d*d); 
         abort();
      }

   }
         
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that fifth powers pass the test */
   {
      bits = n_randint(state, FLINT_BITS/5) + 1;
      d = n_randtest_bits(state, bits);

      result = n_is_perfect_power235(n_pow(d, 5));
      
      if (!result)
      {
         printf("FAIL:\n");
         printf("d^5 = %lu is declared not a perfect power\n", d*d*d*d*d); 
         abort();
      }

   }
         
   for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test that non prefect powers fail */
   {
      mpz_t d_m;
      mpz_init(d_m);

      do
      {
         d = n_randtest(state);
         mpz_set_ui(d_m, d);
      } while (mpz_perfect_power_p(d_m));

      result = !n_is_perfect_power235(d);
      if (!result)
      {
         printf("FAIL:\n");
         printf("d = %lu is declared a perfect power\n", d); 
         abort();
      }

      mpz_clear(d_m);
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
