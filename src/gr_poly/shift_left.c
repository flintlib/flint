/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_shift_left(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i;

    /* Copy in reverse to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = len; i--; )
            status |= gr_set(GR_ENTRY(res, n + i, sz), GR_ENTRY(poly, i, sz), ctx);
    }
    else
    {
        for (i = len; i--; )
            gr_swap(GR_ENTRY(res, n + i, sz), GR_ENTRY(res, i, sz), ctx);
    }

    status |= _gr_vec_zero(res, n, ctx);
    return status;
}

int
gr_poly_shift_left(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (n == 0)
        return gr_poly_set(res, poly, ctx);

    if (poly->length == 0)
        return gr_poly_zero(res, ctx);

    gr_poly_fit_length(res, poly->length + n, ctx);
    status |= _gr_poly_shift_left(res->coeffs, poly->coeffs, poly->length, n, ctx);
    _gr_poly_set_length(res, poly->length + n, ctx);
    return status;
}
