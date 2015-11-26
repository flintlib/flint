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

ulong byte_swap_naive(ulong n)
{
   ulong r = 0;
   slong i;

   for (i = 0; i < sizeof(ulong); i++)
   {
      r <<= 8;
      r |= (n & 0xFF);
      n >>= 8;
   }

   return r;
}

int main(void)
{
   int i, result;
   FLINT_TEST_INIT(state);
   
   flint_printf("byte_swap....");
   fflush(stdout);

   /* byte_swap(byte_swap(n)) == n */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      ulong n, r1, r2;

      n = n_randtest(state);
      
      r1 = n;
      
      r2 = n;
      byte_swap(r2);
      byte_swap(r2);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("byte_swap(byte_swap(n)) != n\n");
         flint_printf("n = %wx, r1 = %wx, r2 = %wx\n", n, r1, r2);
         abort();
      }
   }

   /* byte_swap(n) == byte_swap_naive(n) */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      ulong n, r1, r2;

      n = n_randtest(state);
      
      r1 = n;
      byte_swap(r1);

      r2 = byte_swap_naive(n);

      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("byte_swap(n) != byte_swap_naive(n)\n");
         flint_printf("n = %wx, r1 = %wx, r2 = %wx\n", n, r1, r2);
         abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
