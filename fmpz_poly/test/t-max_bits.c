/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int main(void)
{
   int result;
   printf("max_bits....");
   fflush(stdout);
   
   fmpz_poly_randinit();
   
   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz_poly_t a, b, c;

      fmpz_poly_init(a);
      fmpz_poly_init(b);
      fmpz_poly_init(c);
      ulong bits = n_randint(200);
	  ulong bits2;
	  fmpz_poly_randtest(a, n_randint(100), bits);
      
      bits2 = fmpz_poly_max_bits(a);

      result = (bits >= FLINT_ABS(bits2));
      if (!result)
      {
         printf("Error:\n");
         printf("bits = %ld, bits2 = %ld\n", bits, bits2);
		 abort();
      }

      fmpz_poly_clear(a);
      fmpz_poly_clear(b);
      fmpz_poly_clear(c);
   }

   for (ulong i = 0; i < 10000UL; i++) 
   {
      fmpz_poly_t a;

      fmpz_poly_init(a);
      ulong bits = n_randint(200);
	  ulong bits2;
	  fmpz_poly_randtest_unsigned(a, n_randint(100), bits);
      
      bits2 = fmpz_poly_max_bits(a);

      result = (bits >= bits2);
      if (!result)
      {
         printf("Error:\n");
         printf("bits = %ld, bits2 = %ld\n", bits, bits2);
		 abort();
      }

      fmpz_poly_clear(a);
   }

   fmpz_poly_randclear();
      
   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
