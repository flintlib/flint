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

int _gr_poly_divrem_newton(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    const slong lenQ = lenA - lenB + 1;

    status |= _gr_poly_div_newton(Q, A, lenA, B, lenB, ctx);

    if (lenB > 1 && status == GR_SUCCESS)
    {
        if (R == A)
        {
            gr_ptr W;
            GR_TMP_INIT_VEC(W, lenB - 1, ctx);
            status |= _gr_poly_mullow(W, Q, lenQ, B, lenB - 1, lenB - 1, ctx);
            status |= _gr_vec_sub(R, A, W, lenB - 1, ctx);
            GR_TMP_CLEAR_VEC(W, lenB - 1, ctx);
        }
        else
        {
            status |= _gr_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, ctx);
            status |= _gr_vec_sub(R, A, R, lenB - 1, ctx);
        }
    }

    return status;
}

int
gr_poly_divrem_newton(gr_poly_t Q, gr_poly_t R,
    const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    slong sz = ctx->sizeof_elem;
    gr_poly_t tQ, tR;
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
        gr_poly_init2(tQ, lenQ, ctx);
        q = tQ->coeffs;
    }
    else
    {
        gr_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == B)
    {
        gr_poly_init2(tR, lenB - 1, ctx);
        r = tR->coeffs;
    }
    else
    {
        gr_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    status |= _gr_poly_divrem_newton(q, r, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (Q == A || Q == B)
    {
        gr_poly_swap(tQ, Q, ctx);
        gr_poly_clear(tQ, ctx);
    }
    else
    {
        _gr_poly_set_length(Q, lenQ, ctx);
    }

    if (R == B)
    {
        gr_poly_swap(tR, R, ctx);
        gr_poly_clear(tR, ctx);
    }

    _gr_poly_set_length(Q, lenQ, ctx);
    _gr_poly_set_length(R, lenB - 1, ctx);
    _gr_poly_normalise(R, ctx);

    return status;
}
