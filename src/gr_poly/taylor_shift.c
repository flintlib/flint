/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo: algorithm selection */
int
_gr_poly_taylor_shift_generic(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    if (len <= 20)
        return _gr_poly_taylor_shift_horner(res, poly, len, c, ctx);
    else
        return _gr_poly_taylor_shift_divconquer(res, poly, len, c, ctx);
}

int
gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != f)
        status |= gr_poly_set(res, f, ctx);

    status |= _gr_poly_taylor_shift(res->coeffs, res->coeffs, res->length, c, ctx);
    return status;
}
