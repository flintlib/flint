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
#include "ulong_extras.h"
#include "fmpz.h"

int main(void)
{
   int result;
   printf("bit_pack/bit_unpack....");
   fflush(stdout);

   fmpz_randinit();

   for (ulong i = 0; i < 1000000UL; i++) 
   {
      fmpz_t a, b;
      mp_bitcnt_t bits = n_randint(300) + 1;
      ulong space = (300 - 1)/FLINT_BITS + 2; // 2 to accomodate shift
      mp_limb_t * arr = (mp_limb_t *) calloc(sizeof(mp_limb_t), space);
      mp_bitcnt_t shift = n_randint(FLINT_BITS);
      int negate = (int) n_randint(2);

	  fmpz_init(a);
      fmpz_init(b);
      
      fmpz_randtest(a, bits - 1); // need one bit for sign
      
	  arr[0] = n_randbits(shift);

      fmpz_poly_bit_pack(arr, shift, bits, a, -1, 0);
      fmpz_poly_bit_unpack(b, arr, shift, bits, -1, 0);

      result = (fmpz_cmp(a, b) == 0);

      if (!result)
      {
         printf("FAIL\n");
         fmpz_print(a); printf("\n");
         fmpz_print(b); printf("\n");
		 abort();
      }

      free(arr);
	  fmpz_clear(a);
      fmpz_clear(b);
   }

   fmpz_randclear();

   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
