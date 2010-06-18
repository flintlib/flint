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

ulong fmpz_poly_max_limbs(const fmpz_poly_t poly)
{
	ulong i;
	long limbs, max_limbs = 0;
	
	for (i = 0; i < poly->length; i++)
	{
	   limbs = fmpz_size(poly->coeffs + i);
	   if (limbs > max_limbs)
	      max_limbs = limbs;
	}

	return max_limbs;
}
