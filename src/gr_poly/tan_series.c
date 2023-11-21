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
_gr_poly_tan_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    return _gr_poly_tan_series_newton(res, f, flen, len, 20, ctx);
}

int
gr_poly_tan_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tan_series(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
