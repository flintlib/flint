/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   ulong count = 0;
   mp_limb_t d;
   mpz_t d_m;
   FLINT_TEST_INIT(state);
   
   
   flint_printf("is_prime_pocklington....");
   fflush(stdout);
   
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that primes pass the test */
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest(state) | 1;
         flint_mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = flint_mpz_get_ui(d_m);
      } while (mpz_size(d_m) > 1);

      result = n_is_prime_pocklington(d, 100);
      if (result == -1) count++;

      if (result == 0)
      {
         flint_printf("FAIL:\n");
         flint_printf("d = %wu is declared composite\n", d); 
         abort();
      }

      mpz_clear(d_m);
   }
         
   for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test that composites do not pass */
   {
      mpz_init(d_m);

      do
      {
         d = n_randtest(state) | 1;
         if (d == UWORD(1)) d++;
         flint_mpz_set_ui(d_m, d);
      } while (mpz_probab_prime_p(d_m, 12));

      result = n_is_prime_pocklington(d, 100);
      if (result == -1) count++;

      if (result == 1)
      {
         flint_printf("FAIL:\n");
         flint_printf("d = %wu is declared prime\n", d); 
         abort();
      }

      mpz_clear(d_m);
   }

   if (count > 200 * flint_test_multiplier())
   {
      flint_printf("FAIL: Pocklington-Lehmer failed too many times (%wu times)\n", count);
      abort();
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
