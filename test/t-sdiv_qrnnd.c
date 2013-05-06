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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   flint_rand_t state;
   flint_randinit(state);

   printf("sdiv_qrnnd....");
   fflush(stdout);

   for (i = 0; i < 1000000; i++)
   {
      mp_limb_signed_t d, nh, nl, q, r, ph, pl;

      do 
      {
         d = n_randtest_not_zero(state);
         nh = n_randtest(state);
      } while ((FLINT_ABS(nh) >= FLINT_ABS(d)/2) || (nh == LONG_MIN));

      nl = n_randtest(state);

      sdiv_qrnnd(q, r, nh, nl, d);
      smul_ppmm(ph, pl, d, q);
      if (r < 0L) sub_ddmmss(ph, pl, ph, pl, 0UL, -r);
      else add_ssaaaa(ph, pl, ph, pl, 0UL, r);

      result = ((ph == nh) && (pl == nl));

      if (!result)
      {
         printf("FAIL:\n");
         printf("nh = %lu, nl = %lu, d = %lu\n", nh, nl, d); 
         printf("ph = %lu, pl = %lu\n", ph, pl);
         abort();
      }
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
