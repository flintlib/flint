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
_gr_poly_add(gr_ptr res, gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    slong sz;
    int status;
    slong min = FLINT_MIN(len1, len2);

    status = _gr_vec_add(res, poly1, poly2, min, ctx);

    if (len1 > min)
    {
        sz = ctx->sizeof_elem;
        status |= _gr_vec_set(GR_ENTRY(res, min, sz), GR_ENTRY(poly1, min, sz), len1 - min, ctx);
    }

    if (len2 > min)
    {
        sz = ctx->sizeof_elem;
        status |= _gr_vec_set(GR_ENTRY(res, min, sz), GR_ENTRY(poly2, min, sz), len2 - min, ctx);
    }

    return status;
}

int
gr_poly_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    int status;
    slong max = FLINT_MAX(poly1->length, poly2->length);

    gr_poly_fit_length(res, max, ctx);

    status = _gr_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);

    _gr_poly_set_length(res, max, ctx);
    _gr_poly_normalise(res, ctx);

    return status;
}
