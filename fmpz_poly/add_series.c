/*
    Copyright (C) 2008, 2009, 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_add_series(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2, slong n)
{
    slong len1, len2, max = FLINT_MAX(poly1->length, poly2->length);

    if (n < 0)
       n = 0;
 
    max = FLINT_MIN(max, n);
    len1 = FLINT_MIN(poly1->length, max);
    len2 = FLINT_MIN(poly2->length, max);

    fmpz_poly_fit_length(res, max);

    _fmpz_poly_add(res->coeffs, poly1->coeffs, len1, poly2->coeffs,
                   len2);

    _fmpz_poly_set_length(res, max);
    _fmpz_poly_normalise(res);  
}
