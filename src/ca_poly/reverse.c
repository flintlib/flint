/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"
#include "gr_poly.h"

void
_ca_poly_reverse(ca_ptr res, ca_srcptr poly, slong len, slong n, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    GR_MUST_SUCCEED(_gr_poly_reverse(res, poly, len, n, gr_ctx));
}

void
ca_poly_reverse(ca_poly_t res, const ca_poly_t poly, slong n, ca_ctx_t ctx)
{
    slong len = FLINT_MIN(n, poly->length);

    if (len == 0)
    {
        ca_poly_zero(res, ctx);
        return;
    }

    ca_poly_fit_length(res, n, ctx);
    _ca_poly_reverse(res->coeffs, poly->coeffs, len, n, ctx);
    _ca_poly_set_length(res, n, ctx);
    _ca_poly_normalise(res, ctx);
}
