/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: faster code when possible */
truth_t
gr_poly_is_one(const gr_poly_t poly, gr_ctx_t ctx)
{
    if (poly->length == 0)
    {
        return gr_ctx_is_zero_ring(ctx);
    }
    else if (poly->length == 1)
    {
        return gr_is_one(poly->coeffs, ctx);
    }
    else
    {
        truth_t res1, res2;

        res1 = gr_is_one(poly->coeffs, ctx);
        if (res1 != T_FALSE)
        {
            res2 = _gr_vec_is_zero(gr_poly_coeff_srcptr(poly, 1, ctx), poly->length - 1, ctx);
            res1 = truth_and(res1, res2);
        }

        return res1;
    }
}
