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
_ca_poly_pow_ui(ca_ptr res, ca_srcptr f, slong flen, ulong exp, ca_ctx_t ctx)
{
    _ca_poly_pow_ui_trunc(res, f, flen, exp, exp * (flen - 1) + 1, ctx);
}

void
ca_poly_pow_ui(ca_poly_t res,
    const ca_poly_t poly, ulong exp, ca_ctx_t ctx)
{
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0)
    {
        ca_poly_one(res, ctx);
    }
    else if (flen == 0)
    {
        ca_poly_zero(res, ctx);
    }
    else
    {
        rlen = exp * (flen - 1) + 1;

        if (res != poly)
        {
            ca_poly_fit_length(res, rlen, ctx);
            _ca_poly_pow_ui(res->coeffs, poly->coeffs, flen, exp, ctx);
            _ca_poly_set_length(res, rlen, ctx);
            _ca_poly_normalise(res, ctx);
        }
        else
        {
            ca_poly_t t;
            ca_poly_init2(t, rlen, ctx);
            _ca_poly_pow_ui(t->coeffs, poly->coeffs, flen, exp, ctx);
            _ca_poly_set_length(t, rlen, ctx);
            _ca_poly_normalise(t, ctx);
            ca_poly_swap(res, t, ctx);
            ca_poly_clear(t, ctx);
        }
    }
}
