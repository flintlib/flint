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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_primitive_part(fmpz * rpoly, fmpz_t rden, 
                               const fmpz * poly, const fmpz_t den, len_t len)
{
    _fmpz_poly_primitive_part(rpoly, poly, len);
    fmpz_one(rden);
}

void fmpq_poly_primitive_part(fmpq_poly_t res, const fmpq_poly_t poly)
{
    const len_t len = poly->length;

    if (len == 0)
    {
        fmpq_poly_zero(res);
    }
    else
    {
        fmpq_poly_fit_length(res, len);
        _fmpq_poly_set_length(res, len);
        
        _fmpq_poly_primitive_part(res->coeffs, res->den, 
                                  poly->coeffs, poly->den, len);
    }
}

