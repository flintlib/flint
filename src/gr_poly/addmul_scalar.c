/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "gr_vec.h"

int
gr_poly_addmul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;
    slong sz = ctx->sizeof_elem;

    if (len == 0 || gr_is_zero(c, ctx) == T_TRUE) {
        return GR_SUCCESS;
    }

    status = GR_SUCCESS;

    if (res != poly) {
        gr_poly_fit_length(res, len, ctx);
        if (poly->length > res->length) {
            status |= _gr_vec_zero(GR_ENTRY(res->coeffs, res->length, sz),
                                    poly->length - res->length, ctx);
        }
    }

    status |= _gr_vec_addmul_scalar(res->coeffs, poly->coeffs, len, c, ctx);
    _gr_poly_set_length(res, FLINT_MAX(res->length, poly->length), ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
