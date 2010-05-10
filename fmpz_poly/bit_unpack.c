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

void fmpz_poly_bit_unpack(fmpz_poly_t poly, ulong length, 
						  const mp_limb_t * arr, ulong bit_size, int negate)
{
   mp_bitcnt_t bits = 0;
   mp_size_t limbs = 0;

   ulong l = bit_size/FLINT_BITS;
   ulong b = bit_size%FLINT_BITS;
   ulong i;
   int borrow = 0;

   fmpz_poly_fit_length(poly, length);
   _fmpz_poly_set_length(poly, length);

   for (i = 0; i < length; i++)
   {
      borrow = fmpz_bit_unpack(poly->coeffs + i, arr + limbs, bits, bit_size, negate, borrow);
	  limbs += l;
	  bits += b;
	  if (bits >= FLINT_BITS)
	  {
	     bits -= FLINT_BITS;
		 limbs++;
	  }
   }
}

void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly, ulong length, 
						           const mp_limb_t * arr, ulong bit_size)
{
   mp_bitcnt_t bits = 0;
   mp_size_t limbs = 0;

   ulong l = bit_size/FLINT_BITS;
   ulong b = bit_size%FLINT_BITS;
   ulong i;
   
   fmpz_poly_fit_length(poly, length);
   _fmpz_poly_set_length(poly, length);

   for (i = 0; i < length; i++)
   {
      fmpz_bit_unpack_unsigned(poly->coeffs + i, arr + limbs, bits, bit_size);
	  limbs += l;
	  bits += b;
	  if (bits >= FLINT_BITS)
	  {
	     bits -= FLINT_BITS;
		 limbs++;
	  }
   }
}
