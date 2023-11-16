/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "ca_poly.h"
#include "gr_poly.h"

void
_ca_poly_pow_ui_trunc(ca_ptr res,
    ca_srcptr f, slong flen, ulong exp, slong len, ca_ctx_t ctx)
{
    if (exp <= 2)
    {
        if (exp == 0)
            ca_one(res, ctx);
        else if (exp == 1)
            _ca_vec_set(res, f, len, ctx);
        else
            _ca_poly_mullow(res, f, flen, f, flen, len, ctx);
    }
    else
    {
        gr_ctx_t gr_ctx;
        _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
        GR_MUST_SUCCEED(_gr_poly_pow_series_ui_binexp(res, f, flen, exp, len, gr_ctx));
    }
}

void
ca_poly_pow_ui_trunc(ca_poly_t res,
    const ca_poly_t poly, ulong exp, slong len, ca_ctx_t ctx)
{
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0 && len != 0)
    {
        ca_poly_one(res, ctx);
    }
    else if (flen == 0 || len == 0)
    {
        ca_poly_zero(res, ctx);
    }
    else
    {
        rlen = poly_pow_length(flen, exp, len);

        if (res != poly)
        {
            ca_poly_fit_length(res, rlen, ctx);
            _ca_poly_pow_ui_trunc(res->coeffs,
                poly->coeffs, flen, exp, rlen, ctx);
            _ca_poly_set_length(res, rlen, ctx);
            _ca_poly_normalise(res, ctx);
        }
        else
        {
            ca_poly_t t;
            ca_poly_init2(t, rlen, ctx);
            _ca_poly_pow_ui_trunc(t->coeffs,
                poly->coeffs, flen, exp, rlen, ctx);
            _ca_poly_set_length(t, rlen, ctx);
            _ca_poly_normalise(t, ctx);
            ca_poly_swap(res, t, ctx);
            ca_poly_clear(t, ctx);
        }
    }
}

