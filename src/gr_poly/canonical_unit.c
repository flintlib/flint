/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
gr_poly_canonical_unit(gr_poly_t res,
    const gr_poly_t poly, gr_ctx_t ctx)
{
    if (poly->length == 0)
    {
        return gr_poly_zero(res, ctx);
    }
    else
    {
        gr_srcptr lc;
        gr_ptr c;
        int status;

        lc = gr_poly_entry_srcptr(poly, poly->length - 1, ctx);

        if (gr_is_zero(lc, ctx) == T_FALSE)
        {
            /* todo: avoid the temporary */
            GR_TMP_INIT(c, ctx);
            status = gr_canonical_unit(c, lc, ctx);
            status |= gr_poly_set_scalar(res, c, ctx);
            GR_TMP_CLEAR(c, ctx);
        }
        else
        {
            status = GR_UNABLE;
        }

        return status;
    }
}
