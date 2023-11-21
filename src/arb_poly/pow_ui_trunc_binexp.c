/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

void
_arb_poly_pow_ui_trunc_binexp(arb_ptr res,
    arb_srcptr f, slong flen, ulong exp, slong len, slong prec)
{
    if (exp <= 2)
    {
        if (exp == 0)
            arb_one(res);
        else if (exp == 1)
            _arb_vec_set_round(res, f, len, prec);
        else
            _arb_poly_mullow(res, f, flen, f, flen, len, prec);
    }
    else if (!_arb_vec_is_finite(f, flen))
    {
        _arb_vec_indeterminate(res, len);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_real_arb(ctx, prec);
        GR_MUST_SUCCEED(_gr_poly_pow_series_ui_binexp(res, f, flen, exp, len, ctx));
    }
}

void
arb_poly_pow_ui_trunc_binexp(arb_poly_t res,
    const arb_poly_t poly, ulong exp, slong len, slong prec)
{
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0 && len != 0)
    {
        arb_poly_one(res);
    }
    else if (flen == 0 || len == 0)
    {
        arb_poly_zero(res);
    }
    else
    {
        rlen = poly_pow_length(flen, exp, len);

        if (res != poly)
        {
            arb_poly_fit_length(res, rlen);
            _arb_poly_pow_ui_trunc_binexp(res->coeffs,
                poly->coeffs, flen, exp, rlen, prec);
            _arb_poly_set_length(res, rlen);
            _arb_poly_normalise(res);
        }
        else
        {
            arb_poly_t t;
            arb_poly_init2(t, rlen);
            _arb_poly_pow_ui_trunc_binexp(t->coeffs,
                poly->coeffs, flen, exp, rlen, prec);
            _arb_poly_set_length(t, rlen);
            _arb_poly_normalise(t);
            arb_poly_swap(res, t);
            arb_poly_clear(t);
        }
    }
}

