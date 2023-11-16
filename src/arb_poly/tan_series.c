/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

#define TAN_NEWTON_CUTOFF 20

void
_arb_poly_tan_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);

    hlen = FLINT_MIN(hlen, len);

    if (_gr_poly_tan_series_newton(res, h, hlen, len, TAN_NEWTON_CUTOFF, ctx) != GR_SUCCESS)
        _arb_vec_indeterminate(res, len);
}

void
arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
{
    if (h->length == 0 || n == 0)
    {
        arb_poly_zero(g);
        return;
    }

    arb_poly_fit_length(g, n);
    _arb_poly_tan_series(g->coeffs, h->coeffs, h->length, n, prec);
    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}

