/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
acb_poly_add_series(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong len, slong prec)
{
    slong len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    len1 = FLINT_MIN(len1, len);
    len2 = FLINT_MIN(len2, len);
    len = FLINT_MAX(len1, len2);

    acb_poly_fit_length(res, len);
    _acb_poly_add(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2, prec);
    _acb_poly_set_length(res, len);
    _acb_poly_normalise(res);
}

