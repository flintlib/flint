/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "gr_poly.h"

int
gr_poly_set_fmpz_poly(gr_poly_t res, const fmpz_poly_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, len = src->length;
    gr_ptr res_coeffs;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    gr_poly_fit_length(res, len, ctx);
    res_coeffs = res->coeffs;

    for (i = 0; i < len; i++)
        status |= gr_set_fmpz(GR_ENTRY(res_coeffs, i, sz), src->coeffs + i, ctx);

    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);

    return status;
}
