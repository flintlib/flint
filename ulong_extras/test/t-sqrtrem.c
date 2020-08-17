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
   FLINT_TEST_INIT(state);
   
   flint_printf("sqrtrem....");
   fflush(stdout);

   

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, r1, r2, s1, s2;
      mpz_t a_m, r2_m, s2_m;

      mpz_init(a_m);
      mpz_init(r2_m);
      mpz_init(s2_m);
      
      a = n_randtest(state);
      
      s1 = n_sqrtrem(&r1, a);

      flint_mpz_set_ui(a_m, a);
      mpz_sqrtrem(s2_m, r2_m, a_m);
      r2 = flint_mpz_get_ui(r2_m);
      s2 = flint_mpz_get_ui(s2_m);
      
      result = ((r1 == r2) && (s1 == s2));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, r1 = %wd, r2 = %wu, s1 = %wd, s2 = %wu\n", a, r1, r2, s1, s2); 
         abort();
      }

      mpz_clear(a_m);
      mpz_clear(r2_m);
      mpz_clear(s2_m);
   }

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, r1, r2, s1, s2, bits;
      mpz_t a_m, r2_m, s2_m;

      mpz_init(a_m);
      mpz_init(r2_m);
      mpz_init(s2_m);
      
      bits = n_randint(state, FLINT_BITS/2 + 1);
      a = n_randtest_bits(state, bits);
      a = a*a;
      a += (n_randint(state, 100) - 50);
      s1 = n_sqrtrem(&r1, a);

      flint_mpz_set_ui(a_m, a);
      mpz_sqrtrem(s2_m, r2_m, a_m);
      r2 = flint_mpz_get_ui(r2_m);
      s2 = flint_mpz_get_ui(s2_m);
      
      result = ((r1 == r2) && (s1 == s2));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, r1 = %wd, r2 = %wu, s1 = %wd, s2 = %wu\n", a, r1, r2, s1, s2); 
         abort();
      }

      mpz_clear(a_m);
      mpz_clear(r2_m);
      mpz_clear(s2_m);
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
