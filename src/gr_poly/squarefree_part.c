/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
gr_poly_squarefree_part(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_poly_t t;
    int status = GR_SUCCESS;

    /* todo */
    if (gr_ctx_is_field(ctx) != T_TRUE || gr_ctx_is_finite_characteristic(ctx) != T_FALSE)
        return GR_UNABLE;

    if (poly->length <= 1)
        return gr_poly_one(res, ctx);

    if (poly->length == 2)
    {
        status |= gr_poly_make_monic(res, poly, ctx);
        if (status != GR_SUCCESS)
            return GR_UNABLE;
    }

    gr_poly_init(t, ctx);
    status |= gr_poly_derivative(t, poly, ctx);
    status |= gr_poly_gcd(t, poly, t, ctx);

    if (status == GR_SUCCESS)
    {
        if (t->length == 1)  /* gcd = 1 */
        {
            status |= gr_poly_make_monic(res, poly, ctx);
        }
        else
        {
            /* should be divexact */
            status |= gr_poly_divrem(res, t, poly, t, ctx);
            if (status == GR_SUCCESS)
                status |= gr_poly_make_monic(res, res, ctx);
        }
    }

    gr_poly_clear(t, ctx);

    if (status != GR_SUCCESS)
        status = GR_UNABLE;

    return status;
}
