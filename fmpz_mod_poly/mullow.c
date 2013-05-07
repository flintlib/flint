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

void _fmpz_mod_poly_mullow(fmpz *res, const fmpz *poly1, long len1, 
                                      const fmpz *poly2, long len2, 
                                      const fmpz_t p, long n)
{
    _fmpz_poly_mullow(res, poly1, len1, poly2, len2, n);
    _fmpz_vec_scalar_mod_fmpz(res, res, n, p);
}

void fmpz_mod_poly_mullow(fmpz_mod_poly_t res, 
    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, long n)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;

    if ((len1 == 0) || (len2 == 0) || (n == 0))
    {
        fmpz_mod_poly_zero(res);
        return;
    }

    n = FLINT_MIN(n, len1 + len2 - 1);

    if ((res == poly1) || (res == poly2))
    {
        fmpz *t = _fmpz_vec_init(n);

        if (len1 >= len2)
            _fmpz_mod_poly_mullow(t, poly1->coeffs, len1, 
                                     poly2->coeffs, len2, &(res->p), n);
        else
            _fmpz_mod_poly_mullow(t, poly2->coeffs, len2, 
                                     poly1->coeffs, len1, &(res->p), n);

        _fmpz_vec_clear(res->coeffs, res->alloc);
        res->coeffs = t;
        res->alloc  = n;
        res->length = n;
        _fmpz_mod_poly_normalise(res);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, n);

        if (len1 >= len2)
            _fmpz_mod_poly_mullow(res->coeffs, poly1->coeffs, len1, 
                                               poly2->coeffs, len2, &(res->p), n);
        else
            _fmpz_mod_poly_mullow(res->coeffs, poly2->coeffs, len2, 
                                               poly1->coeffs, len1, &(res->p), n);

        _fmpz_mod_poly_set_length(res, n);
        _fmpz_mod_poly_normalise(res);
    }
}

