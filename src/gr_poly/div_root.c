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
    gr_ptr t, u;
    slong i, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    /* Q_{len-2} = A_{len-1}, Q_{i-1} = A_i + c Q_i, R = A_0 + c Q_0 */

    if (len < 2)
        return gr_zero(R, ctx);

    GR_TMP_INIT2(t, u, ctx);

    if (Q != A)
    {
        status |= gr_set(GR_ENTRY(Q, len - 2, sz), GR_ENTRY(A, len - 1, sz), ctx);

        for (i = len - 2; i > 0; i--)
        {
            status |= gr_mul(GR_ENTRY(Q, i - 1, sz), GR_ENTRY(Q, i, sz), c, ctx);
            status |= gr_add(GR_ENTRY(Q, i - 1, sz), GR_ENTRY(Q, i - 1, sz), GR_ENTRY(A, i, sz), ctx);
        }

        status |= gr_mul(u, Q, c, ctx);
        status |= gr_add(R, A, u, ctx);
    }
    else
    {
        /* In place: t holds the coefficient of A overwritten last, and
           the entries move by swaps (the entry len - 1, beyond Q, is
           left with an unspecified value). */
        gr_swap(t, GR_ENTRY(Q, len - 2, sz), ctx);
        gr_swap(GR_ENTRY(Q, len - 2, sz), GR_ENTRY(Q, len - 1, sz), ctx);

        for (i = len - 2; i > 0; i--)
        {
            /* t = A_i; entry i - 1 holds A_{i-1} */
            gr_swap(t, GR_ENTRY(Q, i - 1, sz), ctx);
            status |= gr_mul(u, GR_ENTRY(Q, i, sz), c, ctx);
            status |= gr_add(GR_ENTRY(Q, i - 1, sz), GR_ENTRY(Q, i - 1, sz), u, ctx);
        }

        status |= gr_mul(u, Q, c, ctx);
        status |= gr_add(t, t, u, ctx);
        gr_swap(R, t, ctx);
    }

    GR_TMP_CLEAR2(t, u, ctx);
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
