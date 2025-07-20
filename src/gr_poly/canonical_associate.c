/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: int _gr_poly_canonical_associate(gr_ptr res, gr_ptr u, gr_srcptr poly, slong len, gr_ctx_t ctx); */

int
gr_poly_canonical_associate(gr_poly_t ux, gr_poly_t u,
    const gr_poly_t poly, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong len = poly->length;

    if (len == 0)
    {
        status |= gr_poly_zero(ux, ctx);
        if (u != NULL)
            status |= gr_poly_one(u, ctx);
    }
    else if (gr_ctx_is_field(ctx) == T_TRUE && u == NULL)
    {
        return gr_poly_make_monic(ux, poly, ctx);
    }
    else
    {
        gr_ptr lc, c;

        if (ux != poly)
            status |= gr_poly_set(ux, poly, ctx);

        FLINT_ASSERT(len == ux->length);

        GR_TMP_INIT(c, ctx);
        lc = gr_poly_coeff_ptr(ux, len - 1, ctx);
        status |= gr_canonical_associate(lc, c, lc, ctx);
        status |= _gr_vec_mul_scalar(ux->coeffs, ux->coeffs, len - 1, c, ctx);
        _gr_poly_normalise(ux, ctx);

        if (u != NULL)
            status |= gr_poly_set_scalar(u, c, ctx);

        GR_TMP_CLEAR(c, ctx);
    }

    return status;
}

