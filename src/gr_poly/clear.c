/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx)
{
    FLINT_ASSERT(poly->length == 0 || gr_is_zero(gr_poly_coeff_srcptr(poly,poly->length - 1, ctx), ctx) != T_TRUE);

    if (poly->coeffs != NULL)
    {
        _gr_vec_clear(poly->coeffs, poly->alloc, ctx);
        flint_free(poly->coeffs);
    }
}
