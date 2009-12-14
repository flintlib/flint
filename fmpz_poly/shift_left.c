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

   Copyright (C) 2008, 2009 William Hart
   
*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_shift_left(fmpz_poly_t res, const fmpz_poly_t poly, const ulong n)
{
   if (n == 0) // special case, no shift
	{
		fmpz_poly_set(res, poly);
		return;
	}
	
	if (poly->length == 0) // nothing to shift
	{
		fmpz_poly_zero(res);
		return;
	}
	
	fmpz_poly_fit_length(res, poly->length + n);
	
	// copy in reverse order to avoid writing over unshifted coeffs
	long i;
	for (i = poly->length - 1; i >= 0; i--) 
		fmpz_set(res->coeffs + i + n, poly->coeffs + i);

   // insert n zeroes
	for (i = 0; i < n; i++) fmpz_zero(res->coeffs + i);
   
   _fmpz_poly_set_length(res, poly->length + n);
}
