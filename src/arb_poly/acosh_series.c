/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

void
_arb_poly_acosh_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        arb_acosh(g, h, prec);
        _arb_vec_zero(g + 1, n - 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_real_arb(ctx, prec);
        if (_gr_poly_acosh_series(g, h, hlen, n, ctx) != GR_SUCCESS)
            _arb_vec_indeterminate(g, n);
    }
}

void
arb_poly_acosh_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
{
    arb_poly_fit_length(g, n);

    /* acosh(0) = i*pi/2 is not real */
    if (h->length == 0 || n == 0)
        _arb_vec_indeterminate(g->coeffs, n);
    else
        _arb_poly_acosh_series(g->coeffs, h->coeffs, h->length, n, prec);

    _arb_poly_set_length(g, n);
    _arb_poly_normalise(g);
}
