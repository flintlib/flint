/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
ca_poly_set_coeff_ca(ca_poly_t poly, slong n, const ca_t x, ca_ctx_t ctx)
{
    ca_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        _ca_vec_zero(poly->coeffs + poly->length, n - poly->length, ctx);
        poly->length = n + 1;
    }

    ca_set(poly->coeffs + n, x, ctx);
    _ca_poly_normalise(poly, ctx);
}
