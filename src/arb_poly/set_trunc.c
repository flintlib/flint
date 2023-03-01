/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"


void
arb_poly_set_trunc(arb_poly_t res, const arb_poly_t poly, slong n)
{
    if (poly == res)
    {
        arb_poly_truncate(res, n);
    }
    else
    {
        slong rlen;

        rlen = FLINT_MIN(n, poly->length);
        while (rlen > 0 && arb_is_zero(poly->coeffs + rlen - 1))
            rlen--;

        arb_poly_fit_length(res, rlen);
        _arb_vec_set(res->coeffs, poly->coeffs, rlen);
        _arb_poly_set_length(res, rlen);
    }
}

