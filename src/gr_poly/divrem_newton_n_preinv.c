/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_divrem_newton_n_preinv(
    gr_ptr Q,
    gr_ptr R,
    gr_srcptr A, slong lenA,
    gr_srcptr B, slong lenB,
    gr_srcptr Binv, slong lenBinv,
    gr_ctx_t ctx)
{
    slong lenQ = lenA - lenB + 1;
    int status = GR_SUCCESS;

    status |= _gr_poly_div_newton_n_preinv(Q, A, lenA, B, lenB, Binv, lenBinv, ctx);

    if (lenB > 1)
    {
        status |= _gr_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, ctx);
        status |= _gr_vec_sub(R, A, R, lenB - 1, ctx);
    }

    return status;
}

int
gr_poly_divrem_newton_n_preinv(gr_poly_t Q,
                                          gr_poly_t R,
                                          const gr_poly_t A,
                                          const gr_poly_t B,
                                          const gr_poly_t Binv,
                                          gr_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenBinv = Binv->length;
    gr_ptr q, r;
    int status = GR_SUCCESS;

    if (lenB == 0)
        return GR_DOMAIN;

    if (lenA < lenB)
    {
        status |= gr_poly_set(R, A, ctx);
        status |= gr_poly_zero(Q, ctx);
        return status;
    }

    if (lenA > 2 * lenB - 2)
        return GR_UNABLE;

    if (Q == A || Q == B || Q == Binv)
    {
        GR_TMP_INIT_VEC(q, lenA - lenB + 1, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, lenA - lenB + 1, ctx);
        q = Q->coeffs;
    }

    if (R == A || R == B || R == Binv)
    {
        GR_TMP_INIT_VEC(r, lenB - 1, ctx);
    }
    else
    {
        gr_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    status |= _gr_poly_divrem_newton_n_preinv(q, r, A->coeffs, lenA,
                                               B->coeffs, lenB, Binv->coeffs,
                                               lenBinv, ctx);

    if (Q == A || Q == B || Q == Binv)
    {
        GR_TMP_CLEAR_VEC(Q->coeffs, Q->alloc, ctx);
        Q->coeffs = q;
        Q->alloc = lenA - lenB + 1;
    }

    if (R == A || R == B || R == Binv)
    {
        GR_TMP_CLEAR_VEC(R->coeffs, R->alloc, ctx);
        R->coeffs = r;
        R->alloc = lenB - 1;
    }

    _gr_poly_set_length(Q, lenA - lenB + 1, ctx);
    _gr_poly_set_length(R, lenB - 1, ctx);
    _gr_poly_normalise(Q, ctx);
    _gr_poly_normalise(R, ctx);

    return status;
}
