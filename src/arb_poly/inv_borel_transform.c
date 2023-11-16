/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, slong len, slong prec)
{
    slong i;

    arb_t t;
    arb_init(t);

    arb_one(t);

    for (i = 0; i < len; i++)
    {
        if (i > 1)
            arb_mul_ui(t, t, i, prec);

        arb_mul(res + i, poly + i, t, prec);
    }

    arb_clear(t);
}

void
arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, slong prec)
{
    arb_poly_fit_length(res, poly->length);
    _arb_poly_inv_borel_transform(res->coeffs, poly->coeffs, poly->length, prec);
    _arb_poly_set_length(res, poly->length);
    _arb_poly_normalise(res);
}

