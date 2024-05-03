/*
    Copyright (C) 2024 Ricardo Buring

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
gr_poly_contains_zero_poly(const gr_poly_t poly, gr_ctx_t ctx)
{
    if (poly->length == 0)
        return T_TRUE;

    return _gr_vec_contains_zero_vec(poly->coeffs, poly->length, ctx);
}
