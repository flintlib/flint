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
   
   flint_printf("ll_mod_preinv....");
   fflush(stdout);

   

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t d, dinv, nh, nl, r1, r2, m;

      d = n_randtest_not_zero(state);
      m = n_randtest(state);
      r1 = n_randtest(state) % d;
      umul_ppmm(nh, nl, m, d);
      add_ssaaaa(nh, nl, nh, nl, UWORD(0), r1);

      dinv = n_preinvert_limb(d);

      r2 = n_ll_mod_preinv(nh, nl, d, dinv);
      
      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("nh = %wu, nl = %wu, d = %wu, dinv = %wu\n", nh, nl, d, dinv); 
         flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
         abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
