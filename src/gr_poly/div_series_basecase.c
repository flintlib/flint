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
_gr_poly_div_series_basecase_noinv(gr_ptr Q,
    gr_srcptr A, slong Alen,
    gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, l;

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
        return status;
    }

    if (len == 2)
    {
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

        return status;
    }

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

    return status;
}

int
_gr_poly_div_series_basecase_preinv1(gr_ptr Q,
    gr_srcptr A, slong Alen,
    gr_srcptr B, slong Blen, gr_srcptr Binv, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, l;
    int is_one;

    if (len == 0)
        return GR_SUCCESS;

    if (Blen == 0)
        return GR_DOMAIN;

    Alen = FLINT_MIN(Alen, len);
    Blen = FLINT_MIN(Blen, len);

    if (Blen == 1)
    {
        status |= _gr_vec_mul_scalar(Q, A, Alen, Binv, ctx);
        status |= _gr_vec_zero(GR_ENTRY(Q, Alen, sz), len - Alen, ctx);
        return status;
    }

    if (len == 2)
    {
        if (Alen == 1)
        {
            status |= gr_mul(Q, A, Binv, ctx);
            status |= gr_mul(GR_ENTRY(Q, 1, sz), Q, Binv, ctx);
            status |= gr_mul(GR_ENTRY(Q, 1, sz), GR_ENTRY(Q, 1, sz), GR_ENTRY(B, 1, sz), ctx);
            status |= gr_neg(GR_ENTRY(Q, 1, sz), GR_ENTRY(Q, 1, sz), ctx);
        }
        else
        {
            status |= gr_mul(Q, A, Binv, ctx);
            status |= gr_mul(GR_ENTRY(Q, 1, sz), Q, GR_ENTRY(B, 1, sz), ctx);
            status |= gr_sub(GR_ENTRY(Q, 1, sz), GR_ENTRY(A, 1, sz), GR_ENTRY(Q, 1, sz), ctx);
            status |= gr_mul(GR_ENTRY(Q, 1, sz), GR_ENTRY(Q, 1, sz), Binv, ctx);
        }

        return status;
    }

    is_one = (gr_is_one(Binv, ctx) == T_TRUE);

    status |= gr_mul(Q, A, Binv, ctx);

    for (i = 1; i < len; i++)
    {
        l = FLINT_MIN(i, Blen - 1);
        status |= _gr_vec_dot_rev(GR_ENTRY(Q, i, sz), (i < Alen) ? GR_ENTRY(A, i, sz) : NULL, 1, GR_ENTRY(B, 1, sz), GR_ENTRY(Q, i - l, sz), l, ctx);

        if (!is_one)
            status |= gr_mul(GR_ENTRY(Q, i, sz), GR_ENTRY(Q, i, sz), Binv, ctx);
    }

    return status;
}

int
_gr_poly_div_series_basecase_generic(gr_ptr Q,
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
        return status;
    }

    if (len == 2)
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

        return status;
    }

    /* useful to detect?  monomial / series */
#if 0
    if (_gr_vec_is_zero(A, Alen - 1, ctx) == T_TRUE)
    {
        status = _gr_poly_inv_series(GR_ENTRY(Q, Alen - 1, sz), B, Blen, len - (Alen - 1), ctx);

        if (status == GR_SUCCESS)
        {
            status |= _gr_vec_mul_scalar(GR_ENTRY(Q, Alen - 1, sz), GR_ENTRY(Q, Alen - 1, sz), len - (Alen - 1), GR_ENTRY(A, Alen - 1, sz), ctx);
            status |= _gr_vec_zero(Q, Alen - 1, ctx);
            return status;
        }
        else
        {
            status = GR_SUCCESS;
        }
    }
#endif

    {
        gr_ptr q;

        GR_TMP_INIT(q, ctx);

        status = gr_inv(q, B, ctx);

        if (status == GR_SUCCESS)
            status = _gr_poly_div_series_basecase_preinv1(Q, A, Alen, B, Blen, q, len, ctx);
        else
            status = _gr_poly_div_series_basecase_noinv(Q, A, Alen, B, Blen, len, ctx);

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
