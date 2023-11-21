/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: public vec function */
void
_gr_vec_reverse_shallow(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);
    slong sz = ctx->sizeof_elem;
    slong i;

    for (i = 0; i < len; i++)
        set_shallow(GR_ENTRY(res, i, sz), GR_ENTRY(vec, len - 1 - i, sz), ctx);
}

int
_gr_poly_div_newton(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong lenQ, lenB2, alloc;
    gr_ptr Arev, Brev;

    lenQ = lenA - lenB + 1;

    alloc = lenQ + FLINT_MIN(lenB, lenQ);
    Arev = GR_TMP_ALLOC(alloc * sz);
    Brev = GR_ENTRY(Arev, lenQ, sz);

    _gr_vec_reverse_shallow(Arev, GR_ENTRY(A, lenA - lenQ, sz), lenQ, ctx);

    if (lenB >= lenQ)
    {
        _gr_vec_reverse_shallow(Brev, GR_ENTRY(B, lenB - lenQ, sz), lenQ, ctx);
        lenB2 = lenQ;
    }
    else
    {
        _gr_vec_reverse_shallow(Brev, B, lenB, ctx);
        lenB2 = lenB;
    }

    status |= _gr_poly_div_series(Q, Arev, lenQ, Brev, lenB2, lenQ, ctx);
    status |= _gr_poly_reverse(Q, Q, lenQ, lenQ, ctx);

    GR_TMP_FREE(Arev, alloc * sz);

    return status;
}

int
gr_poly_div_newton(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
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
        status = _gr_poly_div_newton(t->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, Qlen, ctx);
        status = _gr_poly_div_newton(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
    }

    _gr_poly_set_length(Q, Qlen, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
