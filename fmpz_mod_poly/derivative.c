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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_derivative(fmpz *res, const fmpz *poly, long len, 
                               const fmpz_t p)
{
    long j, k = 1;

	for (j = 1; j < len; j++)
	{
        if (k == 0)
            fmpz_zero(res + (j - 1));
        else if (k == 1)
            fmpz_set(res + (j - 1), poly + j);
        else
        {
            fmpz_mul_ui(res + (j - 1), poly + j, k);
            fmpz_mod(res + (j - 1), res + (j - 1), p);
        }

        if (fmpz_equal_ui(p, ++k))
            k = 0;
	}
}

void fmpz_mod_poly_derivative(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly)
{
    const long len = poly->length;

    if (len < 2)
    {
        fmpz_mod_poly_zero(res);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, len - 1);
        _fmpz_mod_poly_derivative(res->coeffs, poly->coeffs, len, &(res->p));
        _fmpz_mod_poly_set_length(res, len - 1);
        _fmpz_mod_poly_normalise(res);
    }
}

