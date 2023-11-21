/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_inv_borel_transform(acb_ptr res, acb_srcptr poly, slong len, slong prec)
{
    slong i;

    arb_t t;
    arb_init(t);

    arb_one(t);

    for (i = 0; i < len; i++)
    {
        if (i > 1)
            arb_mul_ui(t, t, i, prec);

        acb_mul_arb(res + i, poly + i, t, prec);
    }

    arb_clear(t);
}

void
acb_poly_inv_borel_transform(acb_poly_t res, const acb_poly_t poly, slong prec)
{
    acb_poly_fit_length(res, poly->length);
    _acb_poly_inv_borel_transform(res->coeffs, poly->coeffs, poly->length, prec);
    _acb_poly_set_length(res, poly->length);
    _acb_poly_normalise(res);
}

