/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
gr_poly_inv(gr_poly_t res,
    const gr_poly_t poly, gr_ctx_t ctx)
{
    if (poly->length == 0)
    {
        if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
            return gr_poly_zero(res, ctx);
        else
            return GR_DOMAIN;
    }
    else if (poly->length == 1)
    {
        int status;
        gr_poly_fit_length(res, 1, ctx);
        status = gr_inv(res->coeffs, poly->coeffs, ctx);
        _gr_poly_set_length(res, 1, ctx);
        _gr_poly_normalise(res, ctx);
        return status;
    }
    else
    {
        return GR_DOMAIN;
    }
}
