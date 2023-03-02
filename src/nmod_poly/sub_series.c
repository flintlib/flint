/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2008, 2009, 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void nmod_poly_sub_series(nmod_poly_t res, 
            const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)
{
    slong len1, len2, max = FLINT_MAX(poly1->length, poly2->length);

    if (n < 0)
       n = 0;
 
    max = FLINT_MIN(max, n);
    len1 = FLINT_MIN(poly1->length, max);
    len2 = FLINT_MIN(poly2->length, max);

    nmod_poly_fit_length(res, max);

    _nmod_poly_sub(res->coeffs, poly1->coeffs, len1, 
                                    poly2->coeffs, len2, poly1->mod);

    _nmod_poly_set_length(res, max);
    _nmod_poly_normalise(res);
}

