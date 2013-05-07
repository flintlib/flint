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
#include "qsieve.h"

int main(void)
{
   int i;
   flint_rand_t state;
   
   printf("ll_knuth_schroeppel....");
   fflush(stdout);
 
   flint_randinit(state);

   for (i = 0; i < 10000; i++) /* Test random n */
   {
      mp_limb_t hi = 0, lo;
      qs_t qs_inf;
      mp_bitcnt_t bits;
      
      bits = n_randint(state, 2*FLINT_BITS) + 1;
      if (bits > FLINT_BITS)
      {
          lo = n_randlimb(state);
          hi = n_randbits(state, bits - FLINT_BITS);
      } else
          lo = n_randbits(state, bits);
      
      qsieve_ll_init(qs_inf, hi, lo);
      qsieve_ll_knuth_schroeppel(qs_inf);
      qsieve_ll_clear(qs_inf);
   }
   
   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
