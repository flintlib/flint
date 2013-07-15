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
   
   printf("mod_precomp....");
   fflush(stdout);

   flint_randinit(state);

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t bits, d, n, r1, r2;
      double dpre;

      bits = n_randint(state, FLINT_D_BITS) + 1;
      d = n_randtest_bits(state, bits);
      if (bits <= (FLINT_BITS/2)) n = n_randtest(state) % (d*d);
      else n = n_randtest(state);

      /* must have n < 2^(FLINT_BITS - 1) */
      if (FLINT_BIT_COUNT(n) == FLINT_BITS)
         n >>= 1;

      dpre = n_precompute_inverse(d);

      r1 = n_mod_precomp(n, d, dpre);
      r2 = n%d;

      result = (r1 == r2);
      if (!result)
      {
         printf("FAIL:\n");
         printf("n = %lu, d = %lu, dinv = %g\n", n, d, dpre); 
         printf("r1 = %lu, r2 = %lu\n", r1, r2);
         abort();
      }
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
