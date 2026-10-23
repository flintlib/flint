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

truth_t
_ca_poly_check_equal(ca_srcptr poly1, slong len1,
        ca_srcptr poly2, slong len2, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    return _gr_poly_equal(poly1, len1, poly2, len2, gr_ctx);
}

truth_t
ca_poly_check_equal(const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;

    if (len1 >= len2)
        return _ca_poly_check_equal(poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    else
        return _ca_poly_check_equal(poly2->coeffs, len2, poly1->coeffs, len1, ctx);
}
