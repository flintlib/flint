/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_interpolate_exact(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
{
    return _gr_poly_interpolate_exact_newton(res, xs, ys, len, ctx);
}

int
gr_poly_interpolate_exact(gr_poly_t poly,
                                    const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)
{
    int status;
    slong n = xs->length;

    if (n != ys->length)
        return GR_DOMAIN;

    gr_poly_fit_length(poly, n, ctx);
    status = _gr_poly_interpolate_exact(poly->coeffs, xs->entries, ys->entries, n, ctx);
    _gr_poly_set_length(poly, n, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}


int
_gr_poly_interpolate(gr_ptr res, gr_srcptr xs, gr_srcptr ys, slong len, gr_ctx_t ctx)
{
    return _gr_poly_interpolate_newton(res, xs, ys, len, ctx);
}

int
gr_poly_interpolate(gr_poly_t poly,
                                    const gr_vec_t xs, const gr_vec_t ys, gr_ctx_t ctx)
{
    int status;
    slong n = xs->length;

    if (n != ys->length)
        return GR_DOMAIN;

    gr_poly_fit_length(poly, n, ctx);
    status = _gr_poly_interpolate(poly->coeffs, xs->entries, ys->entries, n, ctx);
    _gr_poly_set_length(poly, n, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}

