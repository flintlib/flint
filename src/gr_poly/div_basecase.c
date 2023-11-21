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
_gr_poly_div_basecase_preinv1(gr_ptr Q,
    gr_srcptr A, slong Alen,
    gr_srcptr B, slong Blen, gr_srcptr invB, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong Qlen;
    int is_one;
    slong i, l;

    if (Blen == 1)
        return _gr_vec_mul_scalar(Q, A, Alen, invB, ctx);

    Qlen = Alen - Blen + 1;
    is_one = (gr_is_one(invB, ctx) == T_TRUE);
    status |= gr_mul(GR_ENTRY(Q, Qlen - 1, sz), GR_ENTRY(A, Alen - 1, sz), invB, ctx);

    for (i = 1; i < Qlen; i++)
    {
        l = FLINT_MIN(i, Blen - 1);
        status |= _gr_vec_dot_rev(GR_ENTRY(Q, Qlen - 1 - i, sz), GR_ENTRY(A, Alen - 1 - i, sz), 1,
            GR_ENTRY(B, Blen - 1 - l, sz), GR_ENTRY(Q, Qlen - i, sz), l, ctx);
        if (!is_one)
            status |= gr_mul(GR_ENTRY(Q, Qlen - 1 - i, sz), GR_ENTRY(Q, Qlen - 1 - i, sz), invB, ctx);
    }

    return status;
}

int
_gr_poly_div_basecase_noinv(gr_ptr Q,
    gr_srcptr A, slong Alen,
    gr_srcptr B, slong Blen, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong Qlen;
    slong i, l;

    if (Blen == 1)
        return _gr_vec_div_scalar(Q, A, Alen, B, ctx);

    Qlen = Alen - Blen + 1;
    status = gr_div(GR_ENTRY(Q, Qlen - 1, sz), GR_ENTRY(A, Alen - 1, sz), GR_ENTRY(B, Blen - 1, sz), ctx);

    for (i = 1; status == GR_SUCCESS && i < Qlen; i++)
    {
        l = FLINT_MIN(i, Blen - 1);
        status |= _gr_vec_dot_rev(GR_ENTRY(Q, Qlen - 1 - i, sz), GR_ENTRY(A, Alen - 1 - i, sz), 1,
            GR_ENTRY(B, Blen - 1 - l, sz), GR_ENTRY(Q, Qlen - i, sz), l, ctx);
        status |= gr_div(GR_ENTRY(Q, Qlen - 1 - i, sz), GR_ENTRY(Q, Qlen - 1 - i, sz), GR_ENTRY(B, Blen - 1, sz), ctx);
    }

    return status;
}

int
_gr_poly_div_basecase(gr_ptr Q,
    gr_srcptr A, slong Alen,
    gr_srcptr B, slong Blen, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_ptr invB;

    GR_TMP_INIT(invB, ctx);

    /* todo: we sometimes want to keep dividing, e.g. over RR with small coefficient */
    status = gr_inv(invB, GR_ENTRY(B, Blen - 1, sz), ctx);

    /* constant is a unit; can multiply by inverse */
    if (status == GR_SUCCESS)
        status = _gr_poly_div_basecase_preinv1(Q, A, Alen, B, Blen, invB, ctx);
    else
        status = _gr_poly_div_basecase_noinv(Q, A, Alen, B, Blen, ctx);

    GR_TMP_CLEAR(invB, ctx);

    return status;
}

int
gr_poly_div_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    slong Alen, Blen, Qlen;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    Alen = A->length;
    Blen = B->length;

    if (Blen == 0)
        return GR_DOMAIN;

    if (gr_is_zero(GR_ENTRY(B->coeffs, Blen - 1, sz), ctx) != T_FALSE)
        return GR_UNABLE;

    if (Alen < Blen)
        return gr_poly_zero(Q, ctx);

    Qlen = Alen - Blen + 1;

    if (Q == A || Q == B)
    {
        gr_poly_t t;
        gr_poly_init2(t, Qlen, ctx);
        status = _gr_poly_div_basecase(t->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, Qlen, ctx);
        status = _gr_poly_div_basecase(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
    }

    _gr_poly_set_length(Q, Qlen, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
