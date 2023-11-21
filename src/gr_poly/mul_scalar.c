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
gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0 || gr_is_zero(c, ctx) == T_TRUE)
        return gr_poly_zero(res, ctx);

    if (res != poly)
    {
        gr_poly_fit_length(res, len, ctx);
        _gr_poly_set_length(res, len, ctx);
    }

    status = _gr_vec_mul_scalar(res->coeffs, poly->coeffs, len, c, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
