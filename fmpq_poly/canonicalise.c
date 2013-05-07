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
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_canonicalise(fmpz * poly, fmpz_t den, long len)
{
    if (*den == 1L)
        return;
    
    if (*den == -1L)
    {
        _fmpz_vec_neg(poly, poly, len);
        fmpz_one(den);
    }
    else if (len == 0)
    {
        fmpz_one(den);
    }
    else
    {
        fmpz_t gcd;
        fmpz_init(gcd);
        _fmpz_vec_content(gcd, poly, len);
        if (*gcd != 1L)
            fmpz_gcd(gcd, gcd, den);
        if (fmpz_sgn(den) < 0)
            fmpz_neg(gcd, gcd);
        if (*gcd != 1L)
        {
            _fmpz_vec_scalar_divexact_fmpz(poly, poly, len, gcd);
            fmpz_divexact(den, den, gcd);
        }
        fmpz_clear(gcd);
    }
}

void fmpq_poly_canonicalise(fmpq_poly_t poly)
{
    _fmpq_poly_normalise(poly);
    _fmpq_poly_canonicalise(poly->coeffs, poly->den, poly->length);
}

