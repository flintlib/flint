/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

truth_t
_gr_poly_is_monic(gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    gr_srcptr lead;
    slong sz = ctx->sizeof_elem;
    truth_t is_one, is_zero;

    if (len == 0)
        return T_FALSE;

    lead = GR_ENTRY(poly, len - 1, sz);

    is_one = gr_is_one(lead, ctx);

    if (is_one == T_TRUE)
        return T_TRUE;

    /* in case of unnormalized input */
    is_zero = gr_is_zero(lead, ctx);
    if (is_one == T_FALSE && is_zero == T_FALSE)
        return T_FALSE;

    return T_UNKNOWN;
}

truth_t
gr_poly_is_monic(const gr_poly_t res, gr_ctx_t ctx)
{
    return _gr_poly_is_monic(res->coeffs, res->length, ctx);
}
