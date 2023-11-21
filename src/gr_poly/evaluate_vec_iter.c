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
_gr_poly_evaluate_vec_iter(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx)
{
    slong i;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    for (i = 0; i < n; i++)
        status |= _gr_poly_evaluate(GR_ENTRY(ys, i, sz), poly, plen, GR_ENTRY(xs, i, sz), ctx);

    return status;
}

int
gr_poly_evaluate_vec_iter(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx)
{
    gr_vec_set_length(ys, xs->length, ctx);
    return _gr_poly_evaluate_vec_iter(ys->entries, poly->coeffs, poly->length, xs->entries, xs->length, ctx);
}
