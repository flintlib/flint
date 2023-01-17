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
_gr_poly_div_series_basecase(gr_ptr Q,
    gr_srcptr A, slong Alen,
    gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len == 0)
        return GR_SUCCESS;

    if (Blen == 0)
        return GR_DOMAIN;

    Alen = FLINT_MIN(Alen, len);
    Blen = FLINT_MIN(Blen, len);

    if (Blen == 1)
    {
        status |= _gr_vec_div_scalar(Q, A, Alen, B, ctx);
        status |= _gr_vec_zero(GR_ENTRY(Q, Alen, sz), len - Alen, ctx);
    }
    else if (len == 2)
    {
        /* todo: in appropriate cases, don't do a division */
        if (Alen == 1)
        {
            status |= gr_div(Q, A, B, ctx);
            status |= gr_div(GR_ENTRY(Q, 1, sz), Q, B, ctx);
            status |= gr_mul(GR_ENTRY(Q, 1, sz), GR_ENTRY(Q, 1, sz), GR_ENTRY(B, 1, sz), ctx);
            status |= gr_neg(GR_ENTRY(Q, 1, sz), GR_ENTRY(Q, 1, sz), ctx);
        }
        else
        {
            status |= gr_div(Q, A, B, ctx);
            status |= gr_mul(GR_ENTRY(Q, 1, sz), Q, GR_ENTRY(B, 1, sz), ctx);
            status |= gr_sub(GR_ENTRY(Q, 1, sz), GR_ENTRY(A, 1, sz), GR_ENTRY(Q, 1, sz), ctx);
            status |= gr_div(GR_ENTRY(Q, 1, sz), GR_ENTRY(Q, 1, sz), B, ctx);
        }
    }
    else
    {
        gr_ptr q;
        int is_one;
        slong i, l;

        GR_TMP_INIT(q, ctx);

        status = gr_inv(q, B, ctx);

        /* constant is a unit; can multiply by inverse */
        if (status == GR_SUCCESS)
        {
            is_one = (gr_is_one(q, ctx) == T_TRUE);

            status |= gr_mul(Q, A, q, ctx);

            for (i = 1; i < len; i++)
            {
                l = FLINT_MIN(i, Blen - 1);
                status |= _gr_vec_dot_rev(GR_ENTRY(Q, i, sz), (i < Alen) ? GR_ENTRY(A, i, sz) : NULL, 1, GR_ENTRY(B, 1, sz), GR_ENTRY(Q, i - l, sz), l, ctx);

                if (!is_one)
                    status |= gr_mul(GR_ENTRY(Q, i, sz), GR_ENTRY(Q, i, sz), q, ctx);
            }
        }
        else  /* need to check divisions */
        {
            status = gr_div(Q, A, B, ctx);

            if (status == GR_SUCCESS)
            {
                for (i = 1; i < len; i++)
                {
                    l = FLINT_MIN(i, Blen - 1);

                    status |= _gr_vec_dot_rev(GR_ENTRY(Q, i, sz), (i < Alen) ? GR_ENTRY(A, i, sz) : NULL, 1, GR_ENTRY(B, 1, sz), GR_ENTRY(Q, i - l, sz), l, ctx);
                    status |= gr_div(GR_ENTRY(Q, i, sz), GR_ENTRY(Q, i, sz), B, ctx);

                    if (status != GR_SUCCESS)
                        break;
                }
            }
        }

        GR_TMP_CLEAR(q, ctx);
    }

    return status;
}

int
gr_poly_div_series_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len == 0)
        return gr_poly_zero(Q, ctx);

    if (B->length == 0)
        return GR_DOMAIN;

    if (A->length == 0)
    {
        truth_t is_zero = gr_poly_is_zero(B, ctx);

        if (is_zero == T_FALSE)
            return gr_poly_zero(Q, ctx);

        return GR_UNABLE;
    }

    if (Q == A || Q == B)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_div_series_basecase(t, A, B, len, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Q, len, ctx);
    status = _gr_poly_div_series_basecase(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, len, ctx);
    _gr_poly_set_length(Q, len, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
