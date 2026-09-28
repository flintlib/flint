/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "gr_poly.h"

void
_acb_poly_integral(acb_ptr res, acb_srcptr poly, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    GR_MUST_SUCCEED(_gr_poly_integral(res, poly, len, ctx));
}

void
acb_poly_integral(acb_poly_t res, const acb_poly_t poly, slong prec)
{
    acb_poly_fit_length(res, poly->length + 1);
    _acb_poly_integral(res->coeffs, poly->coeffs, poly->length + 1, prec);
    _acb_poly_set_length(res, poly->length + 1);
    _acb_poly_normalise(res);
}
