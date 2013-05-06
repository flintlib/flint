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
/*****************************************************************************

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int result = 1;
   ulong i, j;
   mp_limb_t p;
   mpz_t i_m;
   
   printf("compute_primes....");
   fflush(stdout);
   
   n_compute_primes(500UL);
   
   if (flint_num_primes < (500UL))
   {
      printf("FAIL:\n");
      printf("Not enough primes computed, flint_num_primes = %lu\n", flint_num_primes);
      abort();
   }
      
   n_compute_primes(10000UL);
   
   if (flint_num_primes < (10000UL))
   {
      printf("FAIL:\n");
      printf("Not enough primes computed, flint_num_primes = %lu\n", flint_num_primes);
      abort();
   }
      
   n_compute_primes(20000UL);
   
   if (flint_num_primes < (20000UL))
   {
      printf("FAIL:\n");
      printf("Not enough primes computed, flint_num_primes = %lu\n", flint_num_primes);
      abort();
   }
      
   n_compute_primes(74000UL);
   
   if (flint_num_primes < (74000UL))
   {
      printf("FAIL:\n");
      printf("Not enough primes computed, flint_num_primes = %lu\n", flint_num_primes);
      abort();
   }
      
   p = flint_primes[0];
   j = 0;
   mpz_init(i_m);
   
   for (i = 0; i < flint_primes_cutoff; i++) /* Test that primes pass the test */
   {
      mpz_set_ui(i_m, i);
      if (mpz_probab_prime_p(i_m, 20))
      {
         result = (p == i);
         
         if (!result)
         {
            printf("FAIL:\n");
            printf("%lu, %lu\n", i, p); 
            abort();
         }

         j++;
         if (j < flint_num_primes) 
            p = flint_primes[j];
      }
   }

   mpz_clear(i_m);
   printf("PASS\n");
   return 0;
}
