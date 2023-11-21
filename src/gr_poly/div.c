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

/* todo: tuning */
#define DIVCONQUER_CUTOFF 10

int
_gr_poly_div_generic(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    int status;

    if (lenB <= DIVCONQUER_CUTOFF || lenA - lenB <= DIVCONQUER_CUTOFF)
    {
        status = _gr_poly_div_basecase(Q, A, lenA, B, lenB, ctx);
    }
    else
    {
        status = _gr_poly_div_newton(Q, A, lenA, B, lenB, ctx);

        /* Newton requires invertible lc(B); basecase and divide-and-conquer
           may yet succeed without it. */
        if (status == GR_DOMAIN)
            status = _gr_poly_div_divconquer_noinv(Q, A, lenA, B, lenB, DIVCONQUER_CUTOFF, ctx);
    }

    return status;
}

int
gr_poly_div(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
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
        status = _gr_poly_div(t->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, Qlen, ctx);
        status = _gr_poly_div(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
    }

    _gr_poly_set_length(Q, Qlen, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
