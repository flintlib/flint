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

    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev

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
   
   flint_printf("mulmod_shoup....");
   fflush(stdout);

   

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, b, d, r1, r2, q, p1, p2, w_pr;
      
      d = n_randtest_not_zero(state) / 2 + 1;
      a = n_randtest(state) % d;
      b = n_randtest(state) % d;
      
      w_pr = n_mulmod_precomp_shoup(a, d);

      r1 = n_mulmod_shoup(a, b, w_pr, d);

      umul_ppmm(p1, p2, a, b);
      p1 %= d;
      udiv_qrnnd(q, r2, p1, p2, d);
      
      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, b = %wu, d = %wu, w_pr = %wu\n", a, b, d, w_pr); 
         flint_printf("q = %wu, r1 = %wu, r2 = %wu\n", q, r1, r2);
         abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}

