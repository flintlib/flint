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
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_mul(fmpz *res, const fmpz *poly1, len_t len1, 
                                   const fmpz *poly2, len_t len2, const fmpz_t p)
{
    _fmpz_poly_mul(res, poly1, len1, poly2, len2);
    _fmpz_vec_scalar_mod_fmpz(res, res, len1 + len2 - 1, p);
}

void fmpz_mod_poly_mul(fmpz_mod_poly_t res, 
                       const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2)
{
    const len_t len1 = poly1->length;
    const len_t len2 = poly2->length;
    const len_t lenr = len1 + len2 - 1;

    if ((len1 == 0) || (len2 == 0))
    {
        fmpz_mod_poly_zero(res);
        return;
    }

    if ((res == poly1) || (res == poly2))
    {
        fmpz *t = _fmpz_vec_init(lenr);

        if (len1 >= len2)
            _fmpz_mod_poly_mul(t, poly1->coeffs, len1, 
                                  poly2->coeffs, len2, &(res->p));
        else
            _fmpz_mod_poly_mul(t, poly2->coeffs, len2, 
                                  poly1->coeffs, len1, &(res->p));

        _fmpz_vec_clear(res->coeffs, res->alloc);
        res->alloc  = lenr;
        res->length = lenr;
        res->coeffs = t;
    }
    else
    {
        fmpz_mod_poly_fit_length(res, lenr);
    
        if (len1 >= len2)
            _fmpz_mod_poly_mul(res->coeffs, poly1->coeffs, len1, 
                                            poly2->coeffs, len2, &(res->p));
        else
            _fmpz_mod_poly_mul(res->coeffs, poly2->coeffs, len2, 
                                            poly1->coeffs, len1, &(res->p));

        _fmpz_mod_poly_set_length(res, lenr);
    }
    _fmpz_mod_poly_normalise(res);
}

