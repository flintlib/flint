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
   flint_randinit(state);

   printf("jacobi....");
   fflush(stdout);

   for (i = 0; i < 10000 * flint_test_multiplier(); i++) 
   {
      mp_limb_t d;
      mpz_t a_m, d_m;
      mp_limb_signed_t a;
      int r1, r2;

      mpz_init(a_m);
      mpz_init(d_m);
      
      a = n_randtest(state);
      d = n_randtest_not_zero(state) | 1L;
      
      r1 = n_jacobi(a, d);

      mpz_set_si(a_m, a);
      mpz_set_ui(d_m, d);
      r2 = mpz_jacobi(a_m, d_m);
      
      result = (r1 == r2);
      if (!result)
      {
         printf("FAIL:\n");
         printf("a = %lu, d = %lu\n", a, d); 
         abort();
      }

      mpz_clear(a_m);
      mpz_clear(d_m);
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
