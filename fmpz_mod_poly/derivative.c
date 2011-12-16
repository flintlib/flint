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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_derivative(fmpz *res, const fmpz *poly, long len, 
                               const fmpz_t p)
{
    _fmpz_poly_derivative(res, poly, len);
    _fmpz_vec_mod_fmpz(res, res, len - 1, p);
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

