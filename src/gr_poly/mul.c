/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
_gr_poly_mul(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    return _gr_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 1, ctx);
}

int
gr_poly_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, len_out, ctx);
        status = _gr_poly_mul(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len_out, ctx);
        status = _gr_poly_mul(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, ctx);
    }

    _gr_poly_set_length(res, len_out, ctx);
    _gr_poly_normalise(res, ctx);  /* needed for non-integral domains */
    return status;
}
