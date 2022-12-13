/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
_gr_poly_divrem_basecase_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong iQ, iR;

    if (R != A)
        status |= _gr_vec_set(R, A, lenA, ctx);

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        if (gr_is_zero(GR_ENTRY(R, iR, sz), ctx) == T_TRUE)
        {
            status |= gr_zero(GR_ENTRY(Q, iQ, sz), ctx);
        }
        else
        {
            status |= gr_mul(GR_ENTRY(Q, iQ, sz), GR_ENTRY(R, iR, sz), invB, ctx);
            status |= _gr_vec_scalar_submul(GR_ENTRY(R, iQ, sz), B, lenB, GR_ENTRY(Q, iQ, sz), ctx);
        }
    }

    return status;
}

int
_gr_poly_divrem_basecase_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong iQ, iR;

    if (R != A)
        status |= _gr_vec_set(R, A, lenA, ctx);

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        if (gr_is_zero(GR_ENTRY(R, iR, sz), ctx) == T_TRUE)
        {
            status |= gr_zero(GR_ENTRY(Q, iQ, sz), ctx);
        }
        else
        {
            status |= gr_div(GR_ENTRY(Q, iQ, sz), GR_ENTRY(R, iR, sz), GR_ENTRY(B, lenB - 1, sz), ctx);
            if (status != GR_SUCCESS)
                break;
            status |= _gr_vec_scalar_submul(GR_ENTRY(R, iQ, sz), B, lenB, GR_ENTRY(Q, iQ, sz), ctx);
        }
    }

    return status;
}

int
_gr_poly_divrem_basecase(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr invB;

    GR_TMP_INIT(invB, ctx);

    status = gr_inv(invB, GR_ENTRY(B, lenB - 1, sz), ctx);

    if (status == GR_SUCCESS)
        status = _gr_poly_divrem_basecase_preinv(Q, R, A, lenA, B, lenB, invB, ctx);
    else
        status = _gr_poly_divrem_basecase_noinv(Q, R, A, lenA, B, lenB, ctx);

    GR_TMP_CLEAR(invB, ctx);

    return status;
}

int
gr_poly_divrem_basecase(gr_poly_t Q, gr_poly_t R,
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

    if (R == B)
    {
        r = flint_malloc(lenA * sz);
        _gr_vec_init(r, lenA, ctx);
    }
    else
    {
        gr_poly_fit_length(R, lenA, ctx);
        r = R->coeffs;
    }

    status |= _gr_poly_divrem_basecase(q, r, A->coeffs, lenA, B->coeffs, lenB, ctx);

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

    if (R == B)
    {
        _gr_vec_clear(R->coeffs, R->alloc, ctx);
        flint_free(R->coeffs);
        R->coeffs = r;
        R->alloc = lenA;
        R->length = lenA;
    }

    _gr_poly_set_length(R, lenB - 1, ctx);
    _gr_poly_normalise(R, ctx);

    return status;
}
