/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_atan_series(ca_ptr res, ca_srcptr f, slong flen, slong len, ca_ctx_t ctx)
{
    ca_t c;

    flen = FLINT_MIN(flen, len);

    if (CA_IS_SPECIAL(f))
    {
        if (ca_is_unknown(f, ctx))
            _ca_vec_unknown(res, len, ctx);
        else
            _ca_vec_undefined(res, len, ctx);
        return;
    }

    ca_init(c, ctx);
    ca_atan(c, f, ctx);

    if (flen == 1)
    {
        _ca_vec_zero(res + 1, len - 1, ctx);
    }
    else
    {
        ca_ptr t, u;
        slong ulen;

        t = _ca_vec_init(len, ctx);
        u = _ca_vec_init(len, ctx);

        /* atan(h(x)) = integral(h'(x)/(1+h(x)^2)) */
        ulen = FLINT_MIN(len, 2 * flen - 1);
        _ca_poly_mullow(u, f, flen, f, flen, ulen, ctx);
        ca_add_ui(u, u, 1, ctx);

        _ca_poly_derivative(t, f, flen, ctx);
        _ca_poly_div_series(res, t, flen - 1, u, ulen, len, ctx);
        _ca_poly_integral(res, res, len, ctx);

        _ca_vec_clear(t, len, ctx);
        _ca_vec_clear(u, len, ctx);
    }

    ca_swap(res, c, ctx);

    if (ca_check_is_number(res, ctx) != T_TRUE)
    {
        if (ca_is_unknown(res, ctx))
            _ca_vec_unknown(res + 1, len - 1, ctx);
        else
            _ca_vec_undefined(res + 1, len - 1, ctx);
        return;
    }

    ca_clear(c, ctx);
}

void
ca_poly_atan_series(ca_poly_t res, const ca_poly_t f, slong len, ca_ctx_t ctx)
{
    slong flen = f->length;

    if (flen == 0 || len == 0)
    {
        ca_poly_zero(res, ctx);
        return;
    }

    ca_poly_fit_length(res, len, ctx);
    _ca_poly_atan_series(res->coeffs, f->coeffs, flen, len, ctx);
    _ca_poly_set_length(res, len, ctx);
    _ca_poly_normalise(res, ctx);
}
