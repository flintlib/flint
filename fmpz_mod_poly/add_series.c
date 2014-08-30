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
    Copyright (C) 2008, 2009, 2014 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_add_series(fmpz_mod_poly_t res, 
            const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n)
{
    slong len1, len2, max = FLINT_MAX(poly1->length, poly2->length);

    if (n < 0)
       n = 0;
 
    max = FLINT_MIN(max, n);
    len1 = FLINT_MIN(poly1->length, max);
    len2 = FLINT_MIN(poly2->length, max);

    fmpz_mod_poly_fit_length(res, max);

    _fmpz_mod_poly_add(res->coeffs, poly1->coeffs, len1, 
                                    poly2->coeffs, len2, &(res->p));

    _fmpz_mod_poly_set_length(res, max);
    _fmpz_mod_poly_normalise(res);
}

