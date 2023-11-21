/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"


void
_ca_poly_log_series(ca_ptr res, ca_srcptr f, slong flen, slong len, ca_ctx_t ctx)
{
    flen = FLINT_MIN(flen, len);

    if (CA_IS_SPECIAL(f))
    {
        if (ca_is_unknown(f, ctx))
            _ca_vec_unknown(res, len, ctx);
        else
            _ca_vec_undefined(res, len, ctx);
        return;
    }

    if (flen == 1)
    {
        ca_log(res, f, ctx);
        _ca_vec_zero(res + 1, len - 1, ctx);
    }
    else if (len == 2)
    {
        ca_div(res + 1, f + 1, f + 0, ctx);  /* safe since hlen >= 2 */
        ca_log(res, f, ctx);
    }
    else if (_ca_vec_check_is_zero(f + 1, flen - 2, ctx) == T_TRUE)  /* f = a + bx^d */
    {
        slong i, j, d = flen - 1;

        for (i = 1, j = d; j < len; j += d, i++)
        {
            if (i == 1)
                ca_div(res + j, f + d, f + 0, ctx);
            else
                ca_mul(res + j, res + j - d, res + d, ctx);
            _ca_vec_zero(res + j - d + 1, flen - 2, ctx);
        }
        _ca_vec_zero(res + j - d + 1, len - (j - d + 1), ctx);

        for (i = 2, j = 2 * d; j < len; j += d, i++)
            ca_div_si(res + j, res + j, i % 2 ? i : -i, ctx);

        ca_log(res, f, ctx); /* done last to allow aliasing */
    }
    else
    {
        ca_ptr f_diff, f_inv;
        ca_t a;
        slong alloc;

        alloc = len + flen - 1;
        f_inv = _ca_vec_init(alloc, ctx);
        f_diff = f_inv + len;

        ca_init(a, ctx);
        ca_log(a, f, ctx);

        _ca_poly_derivative(f_diff, f, flen, ctx);
        _ca_poly_inv_series(f_inv, f, flen, len, ctx);
        _ca_poly_mullow(res, f_inv, len - 1, f_diff, flen - 1, len - 1, ctx);
        _ca_poly_integral(res, res, len, ctx);
        ca_swap(res, a, ctx);

        ca_clear(a, ctx);
        _ca_vec_clear(f_inv, alloc, ctx);
    }

    if (ca_check_is_number(res, ctx) != T_TRUE)
    {
        if (ca_is_unknown(res, ctx))
            _ca_vec_unknown(res + 1, len - 1, ctx);
        else
            _ca_vec_undefined(res + 1, len - 1, ctx);
        return;
    }
}

void
ca_poly_log_series(ca_poly_t res, const ca_poly_t f, slong len, ca_ctx_t ctx)
{
    if (len == 0)
    {
        ca_poly_zero(res, ctx);
        return;
    }

    ca_poly_fit_length(res, len, ctx);
    if (f->length == 0)
    {
        ca_neg_inf(res->coeffs, ctx);
        _ca_vec_undefined(res->coeffs + 1, len - 1, ctx);
    }
    else
    {
        _ca_poly_log_series(res->coeffs, f->coeffs, f->length, len, ctx);
    }
    _ca_poly_set_length(res, len, ctx);
    _ca_poly_normalise(res, ctx);
}
