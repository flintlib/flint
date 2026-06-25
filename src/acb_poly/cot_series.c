/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "gr_poly.h"

void
_acb_poly_cot_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec)
{
    hlen = FLINT_MIN(hlen, n);

    if (hlen == 1)
    {
        acb_cot(g, h, prec);
        _acb_vec_zero(g + 1, n - 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_complex_acb(ctx, prec);
        if (_gr_poly_cot_series(g, h, hlen, n, ctx) != GR_SUCCESS)
            _acb_vec_indeterminate(g, n);
    }
}

void
acb_poly_cot_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec)
{
    acb_poly_fit_length(g, n);

    if (h->length == 0 || n == 0)
        _acb_vec_indeterminate(g->coeffs, n);
    else
        _acb_poly_cot_series(g->coeffs, h->coeffs, h->length, n, prec);

    _acb_poly_set_length(g, n);
    _acb_poly_normalise(g);
}
