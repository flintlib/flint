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

int
gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (res != src)
    {
        gr_poly_fit_length(res, src->length, ctx);
        status = _gr_vec_set(res->coeffs, src->coeffs, src->length, ctx);
        _gr_poly_set_length(res, src->length, ctx);
    }

    return status;
}
