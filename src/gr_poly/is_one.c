/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo: faster code when possible */
truth_t
gr_poly_is_one(const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_ptr tmp;
    gr_poly_t one;
    truth_t res;

    GR_TMP_INIT(tmp, ctx);

    if (gr_one(tmp, ctx) != GR_SUCCESS)
    {
        res = T_UNKNOWN;
    }
    else
    {
        one->coeffs = tmp;
        one->length = 1;
        one->alloc = 1;

        res = gr_poly_equal(poly, one, ctx);
    }

    GR_TMP_CLEAR(tmp, ctx);

    return res;
}
