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
   ulong count = 0UL;
   flint_rand_t state;
   
   printf("is_strong_probabprime_precomp....");
   fflush(stdout);

   flint_randinit(state);

   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test that primes pass the test */
   {
      mp_limb_t a, d, norm;
      mpz_t d_m;
      double dpre;
      mp_limb_t bits = n_randint(state, FLINT_D_BITS-1) + 2;

      mpz_init(d_m);

      do
      {
         d = n_randbits(state, bits) | 1;
         mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = mpz_get_ui(d_m);
      } while (FLINT_BIT_COUNT(d) > FLINT_D_BITS);
      if (d == 2UL) d++;
         
      for (j = 0; j < 100; j++)
      {
         do a = n_randint(state, d);
         while (a == 0UL);
      
         dpre = n_precompute_inverse(d);
         count_trailing_zeros(norm, d - 1);
         result = n_is_strong_probabprime_precomp(d, dpre, a, (d - 1)>>norm);

         if (!result)
         {
            printf("FAIL:\n");
            printf("a = %lu, d = %lu\n", a, d); 
            abort();
         }
      }

      mpz_clear(d_m);
   }

   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test that not too many composites pass */
   {
      mp_limb_t a, d, norm;
      mpz_t d_m;
      double dpre;
      mp_limb_t bits = n_randint(state, FLINT_D_BITS-3) + 4;

      mpz_init(d_m);

      do
      {
         d = n_randbits(state, bits) | 1;
         mpz_set_ui(d_m, d);
      } while (mpz_probab_prime_p(d_m, 12));

      for (j = 0; j < 100; j++)
      {
         do a = n_randint(state, d);
         while (a == 0UL);
      
         dpre = n_precompute_inverse(d);
         count_trailing_zeros(norm, d - 1);
         result = !n_is_strong_probabprime_precomp(d, dpre, a, (d - 1)>>norm);

         if (!result) count++;
      }

      mpz_clear(d_m);
   }

   if (count > 220 * flint_test_multiplier()) 
   {
      printf("FAIL:\n");
      printf("count = %lu\n", count);
      abort();
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
