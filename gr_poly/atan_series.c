/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_atan_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    gr_ptr c;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    flen = FLINT_MIN(flen, len);

    if (flen == 0)
        return _gr_vec_zero(res, len, ctx);

    GR_TMP_INIT(c, ctx);

    status |= gr_atan(c, f, ctx);

    if (status == GR_SUCCESS)
    {
        if (flen == 1)
        {
            status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), len - 1, ctx);
        }
        else
        {
            gr_ptr t, u;
            slong ulen;

            ulen = FLINT_MIN(len, 2 * flen - 1);

            GR_TMP_INIT_VEC(t, len + ulen, ctx);
            u = GR_ENTRY(t, len, sz);

            /* atan(h(x)) = integral(h'(x)/(1+h(x)^2)) */
            status |= _gr_poly_mullow(u, f, flen, f, flen, ulen, ctx);
            status |= gr_add_ui(u, u, 1, ctx);

            status |= _gr_poly_derivative(t, f, flen, ctx);
            status |= _gr_poly_div_series(res, t, flen - 1, u, ulen, len, ctx);
            status |= _gr_poly_integral(res, res, len, ctx);

            GR_TMP_CLEAR_VEC(t, len + ulen, ctx);
        }

        gr_swap(res, c, ctx);
    }

    GR_TMP_CLEAR(c, ctx);

    return status;
}

int
gr_poly_atan_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
    {
        return gr_poly_zero(res, ctx);
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_atan_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
