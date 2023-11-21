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
acb_poly_set_trunc_round(acb_poly_t res, const acb_poly_t poly, slong n, slong prec)
{
    if (poly == res)
    {
        acb_poly_truncate(res, n);
        _acb_vec_set_round(res->coeffs, res->coeffs, res->length, prec);
    }
    else
    {
        slong rlen;

        rlen = FLINT_MIN(n, poly->length);
        while (rlen > 0 && acb_is_zero(poly->coeffs + rlen - 1))
            rlen--;

        acb_poly_fit_length(res, rlen);
        _acb_vec_set_round(res->coeffs, poly->coeffs, rlen, prec);
        _acb_poly_set_length(res, rlen);
    }
}

