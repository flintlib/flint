/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_sub) (TEMPLATE(T, struct) * res,
                        const TEMPLATE(T, struct) * poly1, slong len1,
                        const TEMPLATE(T, struct) * poly2, slong len2,
                        const TEMPLATE(T, ctx_t) ctx)
{
    const slong min = FLINT_MIN(len1, len2);
    slong i;

    for (i = 0; i < min; i++)
        TEMPLATE(T, sub) (res + i, poly1 + i, poly2 + i, ctx);

    if (poly1 != res)
        for (i = min; i < len1; i++)
            TEMPLATE(T, set) (res + i, poly1 + i, ctx);

    for (i = min; i < len2; i++)
        TEMPLATE(T, neg) (res + i, poly2 + i, ctx);
}

void
TEMPLATE(T, poly_sub) (TEMPLATE(T, poly_t) res,
                       const TEMPLATE(T, poly_t) poly1,
                       const TEMPLATE(T, poly_t) poly2,
                       const TEMPLATE(T, ctx_t) ctx)
{
    const slong max = FLINT_MAX(poly1->length, poly2->length);

    TEMPLATE(T, poly_fit_length) (res, max, ctx);

    _TEMPLATE(T, poly_sub) (res->coeffs, poly1->coeffs, poly1->length,
                            poly2->coeffs, poly2->length, ctx);

    _TEMPLATE(T, poly_set_length) (res, max, ctx);
    _TEMPLATE(T, poly_normalise) (res, ctx);
}


#endif
