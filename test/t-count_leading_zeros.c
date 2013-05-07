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
 
   printf("count_leading_zeros....");
   fflush(stdout);

   for (i = 0; i < 1000000; i++)
   {
      mp_limb_t n;
      unsigned int count;

      n = n_randtest(state);

      count_leading_zeros(count, n);

      result = ((n == 0UL) || (((long)(n << count) < 0L) && (r_shift(n, FLINT_BITS-count) == 0UL)));
      if (!result)
      {
         printf("FAIL:\n");
         printf("n = %lu, count = %u\n", n, count); 
         abort();
      }
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}
