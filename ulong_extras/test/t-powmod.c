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
   
   flint_printf("powmod....");
   fflush(stdout);

   

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
      } while (n_gcd(d, a) != UWORD(1));
      exp = n_randtest(state);
      
      r1 = n_powmod(a, exp, d);

      flint_mpz_set_ui(a_m, a);
      flint_mpz_set_ui(d_m, d);
      if (exp < WORD(0))
      {
         flint_mpz_powm_ui(r2_m, a_m, -exp, d_m);
         mpz_invert(r2_m, r2_m, d_m);
      } else
         flint_mpz_powm_ui(r2_m, a_m, exp, d_m);      
      r2 = flint_mpz_get_ui(r2_m);
      
      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, exp = %wd, d = %wu\n", a, exp, d); 
         flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
         abort();
      }

      mpz_clear(a_m);
      mpz_clear(d_m);
      mpz_clear(r2_m);
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
