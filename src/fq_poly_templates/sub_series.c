/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_sub_series) (TEMPLATE(T, poly_t) res,
                       const TEMPLATE(T, poly_t) poly1,
                       const TEMPLATE(T, poly_t) poly2,
                       slong n, const TEMPLATE(T, ctx_t) ctx)
{
    slong len1, len2, max = FLINT_MAX(poly1->length, poly2->length);

    if (n < 0)
       n = 0;

    max = FLINT_MIN(max, n);
    len1 = FLINT_MIN(poly1->length, max);
    len2 = FLINT_MIN(poly2->length, max);

    TEMPLATE(T, poly_fit_length) (res, max, ctx);

    _TEMPLATE(T, poly_sub) (res->coeffs, poly1->coeffs, len1,
                            poly2->coeffs, len2, ctx);

    _TEMPLATE(T, poly_set_length) (res, max, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);
}


#endif
