/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_evaluate_rectangular(gr_ptr y, gr_srcptr poly,
    slong len, gr_srcptr x, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len <= 2)
    {
        if (len == 0)
            return gr_zero(y, ctx);

        if (len == 1)
            return gr_set(y, poly, ctx);

        status |= gr_mul(y, x, GR_ENTRY(poly, 1, sz), ctx);
        status |= gr_add(y, y, GR_ENTRY(poly, 0, sz), ctx);
        return status;
    }
    else
    {
        slong i, m, r;
        gr_ptr xs;
        gr_ptr s, t, c;

        m = n_sqrt(len) + 1;
        r = (len + m - 1) / m;

        GR_TMP_INIT_VEC(xs, m + 1, ctx);
        GR_TMP_INIT3(s, t, c, ctx);

        status |= _gr_vec_set_powers(xs, x, m + 1, ctx);
        status |= _gr_vec_dot(y, GR_ENTRY(poly, (r - 1) * m, sz), 0, GR_ENTRY(xs, 1, sz),
            GR_ENTRY(poly, (r - 1) * m + 1, sz), len - (r - 1) * m - 1, ctx);

        for (i = r - 2; i >= 0; i--)
        {
            status |= _gr_vec_dot(s, GR_ENTRY(poly, i * m, sz), 0, GR_ENTRY(xs, 1, sz),
                GR_ENTRY(poly, i * m + 1, sz), m - 1, ctx);
            status |= gr_mul(y, y, GR_ENTRY(xs, m, sz), ctx);
            status |= gr_add(y, y, s, ctx);
        }

        GR_TMP_CLEAR_VEC(xs, m + 1, ctx);
        GR_TMP_CLEAR3(s, t, c, ctx);

        return status;
    }
}

int
gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_rectangular(res, f->coeffs, f->length, x, ctx);
}
