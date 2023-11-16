/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
arb_poly_add_si(arb_poly_t res, const arb_poly_t x, slong y, slong prec)
{
    slong len = x->length;

    if (len == 0)
    {
        arb_poly_set_si(res, y);
    }
    else
    {
        arb_poly_fit_length(res, len);

        arb_add_si(res->coeffs, x->coeffs, y, prec);

        if (res != x)
            _arb_vec_set(res->coeffs + 1, x->coeffs + 1, len - 1);

        _arb_poly_set_length(res, len);
        _arb_poly_normalise(res);
    }
}

