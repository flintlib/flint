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

    Copyright (C) 2009, 2011 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

int main(void)
{
   int i, result;
   flint_rand_t state;
   fmpz_t n, t;
   mp_limb_t fac1, fac2, fac;

   printf("ll_factor....");
   fflush(stdout);
 
   flint_randinit(state);

   fmpz_init(n);
   fmpz_init(t);

   for (i = 0; i < 1000; i++) /* Test random n */
   {
      mp_limb_t hi = 0, lo;
      
      fac1 = n_randprime(state, n_randint(state, FLINT_BITS - 2) + 2, 0);
      do {
         fac2 = n_randprime(state, n_randint(state, FLINT_BITS - 2) + 2, 0);
      } while (fac1 == fac2);

      fmpz_set_ui(n, fac1);
      fmpz_mul_ui(n, n, fac2);
      
      fmpz_set_ui(t, 1);
      fmpz_mul_2exp(t, t, FLINT_BITS);

      fmpz_mod(t, n, t);
      lo = fmpz_get_ui(t);
      
      fmpz_fdiv_q_2exp(t, n, FLINT_BITS);
      hi = fmpz_get_ui(t);

      fac = qsieve_ll_factor(hi, lo);

      /* multiplier may be up to 6 bits, hence the limitation */
      result = ((fac == 0 && fmpz_bits(n) > 2*FLINT_BITS - 6) 
          || fac == fac1 || fac == fac2);
      if (!result)
      {
          printf("FAIL: "); fmpz_print(n); printf(" = %ld * %ld\n", fac1, fac2);
          printf("fac = %ld, bits = %ld\n", fac, fmpz_bits(n));
          abort();
      }
   }
   
   fmpz_clear(t);
   fmpz_clear(n);
   
   flint_randclear(state);
   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
