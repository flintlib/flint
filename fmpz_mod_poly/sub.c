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
    Copyright (C) 2008, 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_sub(fmpz *res, const fmpz *poly1, long len1, 
                                   const fmpz *poly2, long len2, const fmpz_t p)
{
    long i, len = FLINT_MAX(len1, len2);

    _fmpz_poly_sub(res, poly1, len1, poly2, len2);

    for (i = 0; i < len; i++)
    {
        if (fmpz_sgn(res + i) < 0)
            fmpz_add(res + i, res + i, p);
    }
}

void fmpz_mod_poly_sub(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2)
{
    long max = FLINT_MAX(poly1->length, poly2->length);

    fmpz_mod_poly_fit_length(res, max);

    _fmpz_mod_poly_sub(res->coeffs, poly1->coeffs, poly1->length, 
                                    poly2->coeffs, poly2->length, &(res->p));

    _fmpz_mod_poly_set_length(res, max);
    _fmpz_mod_poly_normalise(res);
}

