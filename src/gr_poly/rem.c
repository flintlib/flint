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
_gr_poly_rem(gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    int status;
    gr_ptr Q;
    slong lenQ = lenA - lenB + 1;

    if (lenB == 1)
        return GR_SUCCESS;

    GR_TMP_INIT_VEC(Q, lenQ, ctx);
    status = _gr_poly_divrem(Q, R, A, lenA, B, lenB, ctx);
    GR_TMP_CLEAR_VEC(Q, lenQ, ctx);

    return status;
}

int
gr_poly_rem(gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (lenB == 0)
        return GR_DOMAIN;

    if (gr_is_zero(GR_ENTRY(B->coeffs, lenB - 1, sz), ctx) != T_FALSE)
        return GR_UNABLE;

    if (lenA < lenB)
    {
        status |= gr_poly_set(R, A, ctx);
        return status;
    }

    if (R == B)
    {
        gr_poly_t t;
        gr_poly_init2(t, lenB - 1, ctx);
        status = _gr_poly_rem(t->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
        gr_poly_swap(R, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(R, lenB - 1, ctx);
        status = _gr_poly_rem(R->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
    }

    _gr_poly_set_length(R, lenB - 1, ctx);
    _gr_poly_normalise(R, ctx);

    return status;
}
