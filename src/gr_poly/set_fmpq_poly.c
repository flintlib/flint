/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq_poly.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
gr_poly_set_fmpq_poly(gr_poly_t res, const fmpq_poly_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, len = src->length;
    const fmpz * coeffs = src->coeffs;
    gr_ptr res_coeffs;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    /* todo: len == 1 -> set_fmpq */

    gr_poly_fit_length(res, len, ctx);
    res_coeffs = res->coeffs;

    for (i = 0; i < len; i++)
        status |= gr_set_fmpz(GR_ENTRY(res_coeffs, i, sz), coeffs + i, ctx);

    if (!fmpz_is_one(fmpq_poly_denref(src)))
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);

        status |= gr_set_fmpz(t, fmpq_poly_denref(src), ctx);
        status |= gr_inv(t, t, ctx);

        if (status == GR_SUCCESS)
            status |= _gr_vec_mul_scalar(res_coeffs, res_coeffs, len, t, ctx);

        GR_TMP_CLEAR(t, ctx);
    }

    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);

    return status;
}
