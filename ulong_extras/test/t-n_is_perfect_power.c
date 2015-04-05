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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <flint.h>
#include <ulong_extras.h>

int main()
{
   int i, result, j, k;
   mp_limb_t d, r, t;
   FLINT_TEST_INIT(state);

   flint_printf("n_is_perfect_power....");
   fflush(stdout);

   /* check for random number */
   for (i = 0; i < 100000; i++)
   {
       d = n_randlimb(state);
       result = n_is_perfect_power(&r, d);
       if (result && n_pow(r, result) != d)
       {
         flint_printf("FAIL:\n");
         flint_printf("%wu ** %wu != %wu\n", r, result, d);
         abort();
       }
   }

   /* check for possible power of random numbers, which are not
      perfect powers.
   */

   for (j = 0; j < 100000; j++)
   {
       d = n_randtest_not_zero(state);

       if (d == 1) continue;

       while (d < UWORD_MAX && n_is_perfect_power(&r, d) != 0) d += 1;

       i = FLINT_BIT_COUNT(d);

       for (k = 2; k <= FLINT_BITS / i; k++)
       {
           t = n_pow(d, k);
           result = n_is_perfect_power(&r, t);
           if (result != k)
           {
               flint_printf("FAIL:\n");
               flint_printf("%wu ** %wu != %wu ( = %wu ** %wu)\n", r, result, t, d, k);
               abort();
           }
       }
   }

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}


