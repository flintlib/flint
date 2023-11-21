/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo: tuning */
#define DIVCONQUER_CUTOFF 10

/* todo: distinguish rings / fields */
int
_gr_poly_divrem_generic(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    int status;

    if (lenB == 1 || lenA - lenB == 1)
        return _gr_poly_divrem_basecase(Q, R, A, lenA, B, lenB, ctx);

    status = _gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, ctx);

    /* Newton requires invertible lc(B); basecase and divide-and-conquer
       may yet succeed without it. */
    if (status == GR_DOMAIN)
    {
        if (lenB <= DIVCONQUER_CUTOFF || lenA - lenB <= DIVCONQUER_CUTOFF)
            status = _gr_poly_divrem_basecase_noinv(Q, R, A, lenA, B, lenB, ctx);
        else
            status = _gr_poly_divrem_divconquer_noinv(Q, R, A, lenA, B, lenB, DIVCONQUER_CUTOFF, ctx);
    }

    return status;
}

int
gr_poly_divrem(gr_poly_t Q, gr_poly_t R,
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

    status |= _gr_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (Q == A || Q == B)
    {
        gr_poly_swap(tQ, Q, ctx);
        gr_poly_clear(tQ, ctx);
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
