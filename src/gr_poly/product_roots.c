/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo */
#define _gr_poly_mul_monic _gr_poly_mul

/* todo: parallel version; better temporary management */
int
_gr_poly_product_roots(gr_ptr poly, gr_srcptr xs, slong n, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (n == 0)
    {
        status |= gr_one(poly, ctx);
    }
    else if (n == 1)
    {
        status |= gr_neg(poly, xs, ctx);
        status |= gr_one(GR_ENTRY(poly, 1, sz), ctx);
    }
    else if (n == 2)
    {
        status |= gr_mul(poly, xs, GR_ENTRY(xs, 1, sz), ctx);
        status |= gr_add(GR_ENTRY(poly, 1, sz), xs, GR_ENTRY(xs, 1, sz), ctx);
        status |= gr_neg(GR_ENTRY(poly, 1, sz), GR_ENTRY(poly, 1, sz), ctx);
        status |= gr_one(GR_ENTRY(poly, 2, sz), ctx);
    }
    else if (n == 3)
    {
        status |= gr_mul(GR_ENTRY(poly, 1, sz), xs, GR_ENTRY(xs, 1, sz), ctx);
        status |= gr_mul(poly, GR_ENTRY(poly, 1, sz), GR_ENTRY(xs, 2, sz), ctx);
        status |= gr_neg(poly, poly, ctx);
        status |= gr_add(GR_ENTRY(poly, 2, sz), xs, GR_ENTRY(xs, 1, sz), ctx);
        status |= gr_addmul(GR_ENTRY(poly, 1, sz), GR_ENTRY(poly, 2, sz), GR_ENTRY(xs, 2, sz), ctx);
        status |= gr_add(GR_ENTRY(poly, 2, sz), GR_ENTRY(poly, 2, sz), GR_ENTRY(xs, 2, sz), ctx);
        status |= gr_neg(GR_ENTRY(poly, 2, sz), GR_ENTRY(poly, 2, sz), ctx);
        status |= gr_one(GR_ENTRY(poly, 3, sz), ctx);
    }
    else
    {
        slong m = (n + 1) / 2;
        gr_ptr tmp;

        GR_TMP_INIT_VEC(tmp, n + 2, ctx);

        status |= _gr_poly_product_roots(tmp, xs, m, ctx);
        status |= _gr_poly_product_roots(GR_ENTRY(tmp, m + 1, sz), GR_ENTRY(xs, m, sz), n - m, ctx);
        status |= _gr_poly_mul_monic(poly, tmp, m + 1, GR_ENTRY(tmp, m + 1, sz), n - m + 1, ctx);

        GR_TMP_CLEAR_VEC(tmp, n + 2, ctx);
    }

    return status;
}

int
gr_poly_product_roots(gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)
{
    int status;
    slong n = xs->length;

    gr_poly_fit_length(poly, n + 1, ctx);
    status = _gr_poly_product_roots(poly->coeffs, xs->entries, n, ctx);
    _gr_poly_set_length(poly, n + 1, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}

