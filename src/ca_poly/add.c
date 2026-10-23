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
_ca_poly_add(ca_ptr res, ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    GR_MUST_SUCCEED(_gr_poly_add(res, poly1, len1, poly2, len2, gr_ctx));
}

void
ca_poly_add(ca_poly_t res, const ca_poly_t poly1,
              const ca_poly_t poly2, ca_ctx_t ctx)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    ca_poly_fit_length(res, max, ctx);

    _ca_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length, ctx);

    _ca_poly_set_length(res, max, ctx);
    _ca_poly_normalise(res, ctx);
}
