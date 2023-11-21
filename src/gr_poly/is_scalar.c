/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

truth_t
gr_poly_is_scalar(const gr_poly_t poly, gr_ctx_t ctx)
{
    slong len = poly->length;

    if (len <= 1)
        return T_TRUE;

    return _gr_vec_is_zero(GR_ENTRY(poly->coeffs, 1, ctx->sizeof_elem), len - 1, ctx);
}
