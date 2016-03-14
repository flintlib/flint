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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

#define invert_limb_naive(ninv, n)                    \
   do {                                               \
      mp_limb_t dummy;                                \
      udiv_qrnnd (ninv, dummy, ~(n), ~(WORD(0)), n);  \
   } while (0)

int main(void)
{
   int i, result;
   FLINT_TEST_INIT(state);
   
   flint_printf("invert_limb....");
   fflush(stdout);

   for (i = 0; i < 100000 * flint_test_multiplier(); i++)
   {
      mp_limb_t n, ninv1, ninv2;

      n = n_randtest(state);
      n |= (UWORD(1) << (FLINT_BITS - 1));
      
      invert_limb(ninv1, n);
      invert_limb_naive(ninv2, n);

      result = (ninv1 == ninv2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("n = %wx, ninv1 = %wx, ninv2 = %wx\n", n, ninv1, ninv2); 
         flint_abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
