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
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_compose(fmpz * res, const fmpz * poly1, long len1, 
                               const fmpz * poly2, long len2)
{
    if (len1 == 1)
        fmpz_set(res, poly1);
    else if (len2 == 1)
        _fmpz_poly_evaluate_fmpz(res, poly1, len1, poly2);
    else if (len1 <= 4)
        _fmpz_poly_compose_horner(res, poly1, len1, poly2, len2);
    else
        _fmpz_poly_compose_divconquer(res, poly1, len1, poly2, len2);
}

void
fmpz_poly_compose(fmpz_poly_t res, 
                  const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    long lenr;
    
    if (len1 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }
    if (len1 == 1 || len2 == 0)
    {
        fmpz_poly_set_fmpz(res, poly1->coeffs);
        return;
    }
    
    lenr = (len1 - 1) * (len2 - 1) + 1;
    
    if (res != poly1 && res != poly2)
    {
        fmpz_poly_fit_length(res, lenr);
        _fmpz_poly_compose(res->coeffs, poly1->coeffs, len1, 
                                        poly2->coeffs, len2);
        _fmpz_poly_set_length(res, lenr);
        _fmpz_poly_normalise(res);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, lenr);
        _fmpz_poly_compose(t->coeffs, poly1->coeffs, len1,
                                      poly2->coeffs, len2);
        _fmpz_poly_set_length(t, lenr);
        _fmpz_poly_normalise(t);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
