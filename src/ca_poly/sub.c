/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_sub(ca_ptr res, ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, ca_ctx_t ctx)
{
    slong i, min = FLINT_MIN(len1, len2);

    for (i = 0; i < min; i++)
        ca_sub(res + i, poly1 + i, poly2 + i, ctx);

    for (i = min; i < len1; i++)
        ca_set(res + i, poly1 + i, ctx);

    for (i = min; i < len2; i++)
        ca_neg(res + i, poly2 + i, ctx);
}

void
ca_poly_sub(ca_poly_t res, const ca_poly_t poly1,
              const ca_poly_t poly2, ca_ctx_t ctx)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    ca_poly_fit_length(res, max, ctx);

    _ca_poly_sub(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length, ctx);

    _ca_poly_set_length(res, max, ctx);
    _ca_poly_normalise(res, ctx);
}
