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
   
   printf("powmod....");
   fflush(stdout);

   flint_randinit(state);

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t bits, a, d, r1, r2;
      mpz_t a_m, d_m, r2_m;
      mp_limb_signed_t exp;

      mpz_init(a_m);
      mpz_init(d_m);
      mpz_init(r2_m);
      
      bits = n_randint(state, FLINT_D_BITS) + 1;
      d = n_randtest_bits(state, bits);
      do
      {
         a = n_randtest(state) % d;
      } while (n_gcd(d, a) != 1UL);
      exp = n_randtest(state);
      
      r1 = n_powmod(a, exp, d);

      mpz_set_ui(a_m, a);
      mpz_set_ui(d_m, d);
      if (exp < 0L)
      {
         mpz_powm_ui(r2_m, a_m, -exp, d_m);
         mpz_invert(r2_m, r2_m, d_m);
      } else
         mpz_powm_ui(r2_m, a_m, exp, d_m);      
      r2 = mpz_get_ui(r2_m);
      
      result = (r1 == r2);
      if (!result)
      {
         printf("FAIL:\n");
         printf("a = %lu, exp = %ld, d = %lu\n", a, exp, d); 
         printf("r1 = %lu, r2 = %lu\n", r1, r2);
         abort();
      }

      mpz_clear(a_m);
      mpz_clear(d_m);
      mpz_clear(r2_m);
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
