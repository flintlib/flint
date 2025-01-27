/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
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

void
_gr_vec_reverse_shallow(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx);

int
_gr_poly_div_newton_n_preinv(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr FLINT_UNUSED(B), slong lenB, gr_srcptr Binv, slong lenBinv, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong lenQ;
    gr_ptr Arev;

    lenQ = lenA - lenB + 1;

    if (lenBinv == 0)  /* fallback if the caller failed to invert B */
    {
        status = _gr_vec_zero(Q, lenQ, ctx);
    }
    else
    {
        Arev = GR_TMP_ALLOC(lenQ * sz);
        _gr_vec_reverse_shallow(Arev, GR_ENTRY(A, lenA - lenQ, sz), lenQ, ctx);

        status |= _gr_poly_mullow(Q, Arev, lenQ, Binv, FLINT_MIN(lenQ, lenBinv), lenQ, ctx);
        status |= _gr_poly_reverse(Q, Q, lenQ, lenQ, ctx);

        GR_TMP_FREE(Arev, lenQ * sz);
    }

    return status;
}

int
gr_poly_div_newton_n_preinv(gr_poly_t Q,
                                       const gr_poly_t A,
                                       const gr_poly_t B,
                                       const gr_poly_t Binv,
                                       gr_ctx_t ctx)
{
    slong Alen, Blen, Qlen;
    int status = GR_SUCCESS;
    slong lenBinv = Binv->length;

    Alen = A->length;
    Blen = B->length;

    if (Blen == 0)
        return GR_DOMAIN;

    if (Alen < Blen)
        return gr_poly_zero(Q, ctx);

    if (Alen > 2 * Blen - 2)
        return GR_UNABLE;

    Qlen = Alen - Blen + 1;

    if (Q == A || Q == B || Q == Binv)
    {
        gr_poly_t t;
        gr_poly_init2(t, Qlen, ctx);
        status = _gr_poly_div_newton_n_preinv(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, Binv->coeffs, lenBinv, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, Qlen, ctx);
        status = _gr_poly_div_newton_n_preinv(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, Binv->coeffs, lenBinv, ctx);
    }

    _gr_poly_set_length(Q, Qlen, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
