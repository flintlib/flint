/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_mullow_generic(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    return _gr_poly_mullow_classical(res, poly1, len1, poly2, len2, n, ctx);
}

int
gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong n, gr_ctx_t ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;
    n = FLINT_MIN(n, len_out);

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status = _gr_poly_mullow(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, n, ctx);
        status = _gr_poly_mullow(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

