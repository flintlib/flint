/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
_gr_poly_div_root(gr_ptr Q, gr_ptr R, gr_srcptr A, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    gr_ptr r, t;
    slong i, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len < 2)
        return gr_zero(R, ctx);

    GR_TMP_INIT2(r, t, ctx);

    /* (the extra assignments support aliasing of Q and A) */
    status |= gr_set(t, GR_ENTRY(A, len - 2, sz), ctx);
    status |= gr_set(GR_ENTRY(Q, len - 2, sz), GR_ENTRY(A, len - 1, sz), ctx);
    status |= gr_set(r, GR_ENTRY(Q, len - 2, sz), ctx);

    for (i = len - 2; i > 0; i--)
    {
        status |= gr_mul(r, r, c, ctx);
        status |= gr_add(r, r, t, ctx);
        status |= gr_set(t, GR_ENTRY(A, i - 1, sz), ctx);
        status |= gr_set(GR_ENTRY(Q, i - 1, sz), r, ctx);
    }

    status |= gr_mul(r, r, c, ctx);
    status |= gr_add(R, r, t, ctx);

    GR_TMP_CLEAR2(r, t, ctx);
    return status;
}

int
gr_poly_div_root(gr_poly_t Q, gr_ptr R, const gr_poly_t A, gr_srcptr c, gr_ctx_t ctx)
{
    slong len = A->length;
    int status;

    if (len < 2)
    {
        if (len == 1)
            status = gr_set(R, A->coeffs, ctx);
        else
            status = gr_zero(R, ctx);
        status |= gr_poly_zero(Q, ctx);
        return status;
    }

    gr_poly_fit_length(Q, len - 1, ctx);
    status = _gr_poly_div_root(Q->coeffs, R, A->coeffs, len, c, ctx);
    _gr_poly_set_length(Q, len - 1, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
