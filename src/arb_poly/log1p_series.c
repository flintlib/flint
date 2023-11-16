/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

void
_arb_poly_log1p_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    if (_gr_poly_log1p_series(res, f, flen, n, ctx) != GR_SUCCESS)
        _arb_vec_indeterminate(res, n);
}

void
arb_poly_log1p_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
{
    slong flen = f->length;

    if (flen == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    if (flen == 1 /*&& !arb_contains_si(f->coeffs, -1)*/)
        n = 1;

    arb_poly_fit_length(res, n);
    _arb_poly_log1p_series(res->coeffs, f->coeffs, flen, n, prec);
    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}
