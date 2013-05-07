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

    Copyright (C) 2011 William Hart

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
   flint_randinit(state);

   printf("factor_lehman....");
   fflush(stdout);

   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random numbers */
   {
      mp_limb_t n1, n2;

      do
      {
         n1 = n_randtest_bits(state, n_randint(state, FLINT_BITS) + 1);
      } while (n_is_prime(n1) || (n1 < 2UL)
#if FLINT64 /* cannot compute enough primes */
         || (n1 >= 10000000000000000UL)
#endif
         );

      n2 = n_factor_lehman(n1);
      
      result = ((n1%n2) == 0UL && n1 != n2);
      if (!result)
      {
         printf("FAIL:\n");
         printf("n1 = %lu, n2 = %lu\n", n1, n2); 
         abort();
      }
   }
   
   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random products of two primes */
   {
      mp_limb_t n1, n2, n3, n, limit;

#if FLINT64
      limit = 100000000UL - 100UL;
#else
      limit = 65535UL;
#endif

      n1 = n_randtest(state) % (limit + 1);
      n2 = n_randtest(state) % (limit + 1);
          
      n1 = n_nextprime(n1, 1);
      n2 = n_nextprime(n2, 1);

      /* test a specific bug */
#if FLINT64
      if (i == 0)
      {
           n1 = 72528697;
           n2 = 73339073;
      }
#endif

      n = n1*n2;

      n3 = n_factor_lehman(n);
      
      result = ((n%n3) == 0UL && n != n3);
      if (!result)
      {
         printf("FAIL:\n");
         printf("n = %ld, n3 = %lu\n", n, n3);
         abort();
      }
   }
   
   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
