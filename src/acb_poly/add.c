/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_add(acb_ptr res, acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong prec)
{
    slong i, min = FLINT_MIN(len1, len2);

    for (i = 0; i < min; i++)
        acb_add(res + i, poly1 + i, poly2 + i, prec);

    for (i = min; i < len1; i++)
        acb_set_round(res + i, poly1 + i, prec);

    for (i = min; i < len2; i++)
        acb_set_round(res + i, poly2 + i, prec);
}

void
acb_poly_add(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong prec)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    acb_poly_fit_length(res, max);

    _acb_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length, prec);

    _acb_poly_set_length(res, max);
    _acb_poly_normalise(res);
}

