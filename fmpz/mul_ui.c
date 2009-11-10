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

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_mul_ui(fmpz_t f, const fmpz_t g, const ulong x)
{
	fmpz c2 = *g;
	
	if (x == 0)
	{
		fmpz_zero(f);
		return;
	} else if (!COEFF_IS_MPZ(c2)) // coeff2 is small
	{
		mp_limb_t prod[2];
		mp_limb_t uc2 = FLINT_ABS(c2);
		
		// unsigned limb by limb multiply (assembly for most CPU's)
		umul_ppmm(prod[1], prod[0], uc2, x); 
		if (!prod[1]) // result fits in one limb
		{
			fmpz_set_ui(f, prod[0]);
			if (c2 < 0L) fmpz_neg(f, f);
		} else // result takes two limbs
		{	   
   	   __mpz_struct * mpz_ptr = _fmpz_promote(f);
			// two limbs, least significant first, native endian, no nails, stored in prod
         mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
			if (c2 < 0L) mpz_neg(mpz_ptr, mpz_ptr);	
		}
	} else // coeff2 is large
	{     
		__mpz_struct * mpz_ptr = _fmpz_promote(f); // promote without val as if aliased both are large
      mpz_mul_ui(mpz_ptr, COEFF_TO_PTR(c2), x);		
	}
}
