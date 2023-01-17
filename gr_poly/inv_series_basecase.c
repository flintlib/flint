/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_inv_series_basecase(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len == 0)
        return GR_SUCCESS;

    if (Qlen == 0)
        return GR_DOMAIN;

    Qlen = FLINT_MIN(Qlen, len);

    status = gr_inv(Qinv, Q, ctx);

    if (status != GR_SUCCESS)
        return status;

    if (Qlen == 1)
    {
        status |= _gr_vec_zero(GR_ENTRY(Qinv, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        status |= gr_mul(GR_ENTRY(Qinv, 1, sz), Qinv, Qinv, ctx);
        status |= gr_mul(GR_ENTRY(Qinv, 1, sz), GR_ENTRY(Qinv, 1, sz), GR_ENTRY(Q, 1, sz), ctx);
        status |= gr_neg(GR_ENTRY(Qinv, 1, sz), GR_ENTRY(Qinv, 1, sz), ctx);
    }
    else if ((Qlen == 2 || _gr_vec_is_zero(GR_ENTRY(Q, 1, sz), Qlen - 2, ctx) == T_TRUE))
    {
        /* Special-case binomials */
        /* todo: implement using vector functions (geometric series + inflate) */
        slong i, j, step;

        step = Qlen - 1;

        if (gr_is_one(Qinv, ctx) == T_TRUE)
        {
            status |= gr_neg(GR_ENTRY(Qinv, step, sz), GR_ENTRY(Q, step, sz), ctx);
            for (i = 2 * step; i < len; i += step)
                status |= gr_mul(GR_ENTRY(Qinv, i, sz), GR_ENTRY(Qinv, i - step, sz), GR_ENTRY(Qinv, step, sz), ctx);
        }
        else if (gr_is_neg_one(Qinv, ctx) == T_TRUE)
        {
            status |= gr_neg(GR_ENTRY(Qinv, step, sz), GR_ENTRY(Q, step, sz), ctx);
            for (i = 2 * step; i < len; i += step)
                status |= gr_mul(GR_ENTRY(Qinv, i, sz), GR_ENTRY(Qinv, i - step, sz), GR_ENTRY(Q, step, sz), ctx);
        }
        else
        {
            gr_ptr t;
            GR_TMP_INIT(t, ctx);

            status |= gr_mul(t, Qinv, GR_ENTRY(Q, step, sz), ctx);
            status |= gr_neg(t, t, ctx);

            for (i = step; i < len; i += step)
                status |= gr_mul(GR_ENTRY(Qinv, i, sz), GR_ENTRY(Qinv, i - step, sz), t, ctx);

            GR_TMP_CLEAR(t, ctx);
        }

        for (i = 0; i < len; i += step)
            for (j = i + 1; j < FLINT_MIN(len, i + step); j++)
                status |= gr_zero(GR_ENTRY(Qinv, j, sz), ctx);
    }
    else
    {
        int is_one;
        slong i, l;

        is_one = (gr_is_one(Qinv, ctx) == T_TRUE);

        for (i = 1; i < len; i++)
        {
            l = FLINT_MIN(i, Qlen - 1);

            status |= _gr_vec_dot_rev(GR_ENTRY(Qinv, i, sz), NULL, 1, GR_ENTRY(Q, 1, sz), GR_ENTRY(Qinv, i - l, sz), l, ctx);

            if (!is_one)
                status |= gr_mul(GR_ENTRY(Qinv, i, sz), GR_ENTRY(Qinv, i, sz), Qinv, ctx);
        }
    }

    return status;
}

int
gr_poly_inv_series_basecase(gr_poly_t Qinv, const gr_poly_t Q, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong Qlen;

    if (len == 0)
        return gr_poly_zero(Qinv, ctx);

    Qlen = Q->length;

    if (Qlen == 0)
        return GR_DOMAIN;

    if (Qlen == 1)
        len = 1;

    if (Qinv == Q)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_inv_series_basecase(t, Q, len, ctx);
        gr_poly_swap(Qinv, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Qinv, len, ctx);
    status |= _gr_poly_inv_series_basecase(Qinv->coeffs, Q->coeffs, Q->length, len, ctx);
    _gr_poly_set_length(Qinv, len, ctx);
    _gr_poly_normalise(Qinv, ctx);
    return status;
}
