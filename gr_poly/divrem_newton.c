/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
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

int _gr_poly_divrem_newton(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    const slong lenQ = lenA - lenB + 1;

    status |= _gr_poly_div_newton(Q, A, lenA, B, lenB, ctx);

    if (lenB > 1 && status == GR_SUCCESS)
    {
        status |= _gr_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, ctx);
        status |= _gr_vec_sub(R, A, R, lenB - 1, ctx);
    }

    return status;
}

int
gr_poly_divrem_newton(gr_poly_t Q, gr_poly_t R,
    const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    slong sz = ctx->sizeof_elem;
    gr_ptr q, r;
    int status = GR_SUCCESS;

    if (lenB == 0)
        return GR_DOMAIN;

    if (gr_is_zero(GR_ENTRY(B->coeffs, lenB - 1, sz), ctx) != T_FALSE)
        return GR_UNABLE;

    if (lenA < lenB)
    {
        status |= gr_poly_set(R, A, ctx);
        status |= gr_poly_zero(Q, ctx);
        return status;
    }

    if (Q == A || Q == B)
    {
        q = flint_malloc(lenQ * sz);
        _gr_vec_init(q, lenQ, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = flint_malloc((lenB - 1) * sz);
        _gr_vec_init(r, lenB - 1, ctx);
    }
    else
    {
        gr_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    status |= _gr_poly_divrem_newton(q, r, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (Q == A || Q == B)
    {
        _gr_vec_clear(Q->coeffs, Q->alloc, ctx);
        flint_free(Q->coeffs);
        Q->coeffs = q;
        Q->alloc = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _gr_poly_set_length(Q, lenQ, ctx);
    }

    if (R == A || R == B)
    {
        _gr_vec_clear(R->coeffs, R->alloc, ctx);
        flint_free(R->coeffs);
        R->coeffs = r;
        R->alloc = lenB - 1;
        R->length = lenB - 1;
    }

    _gr_poly_set_length(R, lenB - 1, ctx);
    _gr_poly_normalise(R, ctx);

    return status;
}
