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

   Copyright (C) 2010 William Hart
   
*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void _fmpz_poly_bit_pack(mp_limb_t * arr, const fmpz * poly, ulong len, ulong bit_size, int negate)
{
   mp_bitcnt_t bits = 0;
   mp_size_t limbs = 0;

   ulong l = bit_size/FLINT_BITS;
   ulong b = bit_size%FLINT_BITS;
   ulong i;
   int borrow = 0;

   for (i = 0; i < len; i++)
   {
      borrow = fmpz_bit_pack(arr + limbs, bits, bit_size, poly + i, negate, borrow);
	  limbs += l;
	  bits += b;
	  if (bits >= FLINT_BITS)
	  {
	     bits -= FLINT_BITS;
		 limbs++;
	  }
   }
}
