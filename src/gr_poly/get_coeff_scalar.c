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
gr_poly_get_coeff_scalar(gr_ptr res, const gr_poly_t poly, slong i, gr_ctx_t ctx)
{
    if (i < 0 || i >= poly->length)
        return gr_zero(res, ctx);
    else
        return gr_set(res, GR_ENTRY(poly->coeffs, i, ctx->sizeof_elem), ctx);
}
