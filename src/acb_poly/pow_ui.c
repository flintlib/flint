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
_acb_poly_pow_ui(acb_ptr res, acb_srcptr f, slong flen, ulong exp, slong prec)
{
    _acb_poly_pow_ui_trunc_binexp(res, f, flen, exp, exp * (flen - 1) + 1, prec);
}

void
acb_poly_pow_ui(acb_poly_t res,
    const acb_poly_t poly, ulong exp, slong prec)
{
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0)
    {
        acb_poly_one(res);
    }
    else if (flen == 0)
    {
        acb_poly_zero(res);
    }
    else
    {
        rlen = exp * (flen - 1) + 1;

        if (res != poly)
        {
            acb_poly_fit_length(res, rlen);
            _acb_poly_pow_ui(res->coeffs,
                poly->coeffs, flen, exp, prec);
            _acb_poly_set_length(res, rlen);
            _acb_poly_normalise(res);
        }
        else
        {
            acb_poly_t t;
            acb_poly_init2(t, rlen);
            _acb_poly_pow_ui(t->coeffs,
                poly->coeffs, flen, exp, prec);
            _acb_poly_set_length(t, rlen);
            _acb_poly_normalise(t);
            acb_poly_swap(res, t);
            acb_poly_clear(t);
        }
    }
}

