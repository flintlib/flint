/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "gr_poly.h"

void
_acb_poly_log_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    if (_gr_poly_log_series(res, f, flen, n, ctx) != GR_SUCCESS)
        _acb_vec_indeterminate(res, n);
}

void
acb_poly_log_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, n);
    if (f->length == 0)
        _acb_vec_indeterminate(res->coeffs, n);
    else
        _acb_poly_log_series(res->coeffs, f->coeffs, f->length, n, prec);
    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

