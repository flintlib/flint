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

#define DEFAULT_CUTOFF 8

int
_gr_poly_div_series_generic(gr_ptr Q,
    gr_srcptr A, slong Alen,
    gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx)
{
    int status;

    /* todo */
    if (FLINT_MIN(Blen, len) <= DEFAULT_CUTOFF || ctx->methods[GR_METHOD_POLY_MULLOW] == (gr_funcptr) _gr_poly_mullow_generic)
    {
        status = _gr_poly_div_series_basecase(Q, A, Alen, B, Blen, len, ctx);
    }
    else
    {
        status = _gr_poly_div_series_newton(Q, A, Alen, B, Blen, len, DEFAULT_CUTOFF, ctx);

        /* Newton requires invertible constant term of B; basecase and divide-and-conquer
           may yet succeed without it. */
        if (status == GR_DOMAIN)
            status = _gr_poly_div_series_divconquer(Q, A, Alen, B, Blen, len, DEFAULT_CUTOFF, ctx);
    }

    return status;
}

int
gr_poly_div_series(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
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
        status = gr_poly_div_series(t, A, B, len, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Q, len, ctx);
    status = _gr_poly_div_series(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, len, ctx);
    _gr_poly_set_length(Q, len, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
