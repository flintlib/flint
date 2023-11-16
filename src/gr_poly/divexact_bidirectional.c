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

static int
__gr_poly_divexact_bidirectional(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, int basecase, gr_ctx_t ctx)
{
    slong lenQ, len_lo, len_hi;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    /* See if we can divide from the bottom. */
    while (lenB > 1)
    {
        truth_t is_zero;

        is_zero = gr_is_zero(B, ctx);

        if (is_zero == T_FALSE)
            break;

        /* We cannot tell if the low coefficient of B is zero,
           so series division will not work. However, dividing
           from the top might still be possible. */
        if (is_zero == T_UNKNOWN)
        {
            if (basecase)
                return _gr_poly_divexact_basecase(Q, A, lenA, B, lenB, ctx);
            else
                return _gr_poly_div(Q, A, lenA, B, lenB, ctx);
        }

        /* Discard trailing zero coefficient. B is assumed to divide A,
           there is a corresponding zero coefficient in A,
           which we don't need to check since the contract for divexact
           allows returning nonsense in case of an inexact division. */
        B = GR_ENTRY(B, 1, sz);
        A = GR_ENTRY(A, 1, sz);
        lenB--;
        lenA--;
    }

    if (lenB == 1)
    {
        return _gr_vec_divexact_scalar(Q, A, lenA, B, ctx);
    }

    lenQ = lenA - lenB + 1;
    len_hi = lenQ / 2;
    len_lo = lenQ - len_hi;

    if (basecase)
    {
        status |= _gr_poly_divexact_series_basecase(Q, A, lenA, B, lenB, len_lo, ctx);
        status |= _gr_poly_divexact_basecase(GR_ENTRY(Q, len_lo, sz), GR_ENTRY(A, len_lo, sz), lenA - len_lo, B, lenB, ctx);
    }
    else
    {
        status |= _gr_poly_div_series(Q, A, lenA, B, lenB, len_lo, ctx);
        status |= _gr_poly_div(GR_ENTRY(Q, len_lo, sz), GR_ENTRY(A, len_lo, sz), lenA - len_lo, B, lenB, ctx);
    }

    return status;
}

int
_gr_poly_divexact_bidirectional(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    return __gr_poly_divexact_bidirectional(Q, A, lenA, B, lenB, 0, ctx);
}

int
_gr_poly_divexact_basecase_bidirectional(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    return __gr_poly_divexact_bidirectional(Q, A, lenA, B, lenB, 1, ctx);
}

int
gr_poly_divexact_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
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
        status = _gr_poly_divexact_bidirectional(t->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, Qlen, ctx);
        status = _gr_poly_divexact_bidirectional(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
    }

    _gr_poly_set_length(Q, Qlen, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}

int
gr_poly_divexact_basecase_bidirectional(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
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
        status = _gr_poly_divexact_bidirectional(t->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, Qlen, ctx);
        status = _gr_poly_divexact_bidirectional(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
    }

    _gr_poly_set_length(Q, Qlen, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
