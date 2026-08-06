/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "arb_poly.h"

void
_arb_poly_cot_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        arb_cot_pi(g, h, prec);
        _arb_vec_zero(g + 1, len - 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_real_arb(ctx, prec);
        if (_gr_poly_cot_pi_series(g, h, hlen, len, ctx) != GR_SUCCESS)
            _arb_vec_indeterminate(g, len);
    }
}

void
arb_poly_cot_pi_series(arb_poly_t res, const arb_poly_t f, slong len, slong prec)
{
    arb_poly_fit_length(res, len);

    if (f->length == 0 || len == 0)
        _arb_vec_indeterminate(res->coeffs, len);
    else
        _arb_poly_cot_pi_series(res->coeffs, f->coeffs, f->length, len, prec);

    _arb_poly_set_length(res, len);
    _arb_poly_normalise(res);
}
