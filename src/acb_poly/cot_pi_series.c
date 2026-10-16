/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "acb_poly.h"

void
_acb_poly_cot_pi_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec)
{
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        acb_cot_pi(g, h, prec);
        _acb_vec_zero(g + 1, len - 1);
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_complex_acb(ctx, prec);
        if (_gr_poly_cot_pi_series(g, h, hlen, len, ctx) != GR_SUCCESS)
            _acb_vec_indeterminate(g, len);
    }
}

void
acb_poly_cot_pi_series(acb_poly_t res, const acb_poly_t f, slong len, slong prec)
{
    acb_poly_fit_length(res, len);

    if (f->length == 0 || len == 0)
        _acb_vec_indeterminate(res->coeffs, len);
    else
        _acb_poly_cot_pi_series(res->coeffs, f->coeffs, f->length, len, prec);

    _acb_poly_set_length(res, len);
    _acb_poly_normalise(res);
}
