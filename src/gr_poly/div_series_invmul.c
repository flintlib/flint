/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2019 William Hart
    Copyright (C) 2014, 2021, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_div_series_invmul(gr_ptr res, gr_srcptr B, slong Blen, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx)
{
    gr_ptr Ainv;
    int status = GR_SUCCESS;

    GR_TMP_INIT_VEC(Ainv, len, ctx);

    status |= _gr_poly_inv_series(Ainv, A, Alen, len, ctx);
    if (status == GR_SUCCESS)
        status |= _gr_poly_mullow(res, Ainv, len, B, Blen, len, ctx);

    GR_TMP_CLEAR_VEC(Ainv, len, ctx);

    return status;
}

int
gr_poly_div_series_invmul(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx)
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
        status = gr_poly_div_series_invmul(t, A, B, len, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Q, len, ctx);
    status = _gr_poly_div_series_invmul(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, len, ctx);
    _gr_poly_set_length(Q, len, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
