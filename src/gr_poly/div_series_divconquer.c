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
_gr_poly_div_series_divconquer(gr_ptr res, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, slong cutoff, gr_ctx_t ctx)
{
    gr_ptr Arev, Brev;
    slong Arevlen;
    int status = GR_SUCCESS;

    Alen = FLINT_MIN(Alen, len);
    Blen = FLINT_MIN(Blen, len);

    Arevlen = Blen + len - 1;

    /* todo: algorithm without zero padding */
    /* todo: shallow reversals */
    GR_TMP_INIT_VEC(Arev, Arevlen, ctx);
    GR_TMP_INIT_VEC(Brev, Blen, ctx);

    status |= _gr_poly_reverse(Arev, A, Alen, Arevlen, ctx);
    status |= _gr_poly_reverse(Brev, B, Blen, Blen, ctx);
    status |= _gr_poly_div_divconquer(res, Arev, Arevlen, Brev, Blen, cutoff, ctx);
    status |= _gr_poly_reverse(res, res, len, len, ctx);

    GR_TMP_CLEAR_VEC(Arev, Arevlen, ctx);
    GR_TMP_CLEAR_VEC(Brev, Blen, ctx);

    return status;
}

int
gr_poly_div_series_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)
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
        status = gr_poly_div_series_divconquer(t, A, B, len, cutoff, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Q, len, ctx);
    status = _gr_poly_div_series_divconquer(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, len, cutoff, ctx);
    _gr_poly_set_length(Q, len, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
