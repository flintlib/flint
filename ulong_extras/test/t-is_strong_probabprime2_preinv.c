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
   ulong count = UWORD(0);
   FLINT_TEST_INIT(state);
   
   flint_printf("is_strong_probabprime2_preinv....");
   fflush(stdout);

   

   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test that primes pass the test */
   {
      mp_limb_t a, d, dinv, norm;
      mpz_t d_m;
      ulong j;
      
      mpz_init(d_m);

      do
      {
         d = n_randtest(state) | 1;
         flint_mpz_set_ui(d_m, d);
         mpz_nextprime(d_m, d_m);
         d = flint_mpz_get_ui(d_m);
      } while (mpz_size(d_m) > 1);
      if (d == UWORD(2)) d++;
         
      for (j = 0; j < 100; j++)
      {
         do a = n_randtest(state) % d;
         while (a == UWORD(0));
      
         dinv = n_preinvert_limb(d);
         count_trailing_zeros(norm, d - 1);
         result = n_is_strong_probabprime2_preinv(d, dinv, a, (d - 1)>>norm);

         if (!result)
         {
            flint_printf("FAIL:\n");
            flint_printf("a = %wu, d = %wu\n", a, d); 
            fflush(stdout);
            flint_abort();
         }
      }

      mpz_clear(d_m);
   }
         
   for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test that not too many composites pass */
   {
      mp_limb_t a, d, dinv, norm;
      mpz_t d_m;
      ulong j;
      
      mpz_init(d_m);

      do
      {
         d = n_randtest(state) | 1;
         if (d == UWORD(1)) d++;
         flint_mpz_set_ui(d_m, d);
      } while (mpz_probab_prime_p(d_m, 12));

      for (j = 0; j < 100; j++)
      {
         do a = n_randtest(state) % d;
         while (a == UWORD(0));
      
         dinv = n_preinvert_limb(d);
         count_trailing_zeros(norm, d - 1);
         result = !n_is_strong_probabprime2_preinv(d, dinv, a, (d - 1)>>norm);

         if (!result) count++;
      }

      mpz_clear(d_m);
   }

#if FLINT64
   if (count > 220 * flint_test_multiplier()) 
#else
   if (count > 432 * flint_test_multiplier())
#endif
   {
      flint_printf("FAIL:\n");
      flint_printf("count = %wu\n", count);
      fflush(stdout);
      flint_abort();
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
