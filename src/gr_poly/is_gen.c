/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo: faster code when possible */
truth_t
gr_poly_is_gen(const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_ptr tmp;
    gr_poly_t gen;
    truth_t res;
    slong sz = ctx->sizeof_elem;

    GR_TMP_INIT_VEC(tmp, 2, ctx);

    if (gr_one(GR_ENTRY(tmp, 1, sz), ctx) != GR_SUCCESS)
    {
        res = T_UNKNOWN;
    }
    else
    {
        res = gr_is_zero(GR_ENTRY(tmp, 1, sz), ctx);

        if (res != T_UNKNOWN)
        {
            gen->coeffs = tmp;
            gen->length = (res == T_TRUE) ? 1 : 2;
            gen->alloc = gen->length;
            res = gr_poly_equal(poly, gen, ctx);
        }
    }

    GR_TMP_CLEAR_VEC(tmp, 2, ctx);

    return res;
}
