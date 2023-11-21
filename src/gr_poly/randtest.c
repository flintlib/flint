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
gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_poly_fit_length(poly, len, ctx);
    status |= _gr_vec_randtest(poly->coeffs, state, len, ctx);
    _gr_poly_set_length(poly, len, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}
