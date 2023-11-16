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
ca_poly_set_ca(ca_poly_t poly, const ca_t x, ca_ctx_t ctx)
{
    if (ca_check_is_zero(x, ctx) == T_TRUE)
    {
        ca_poly_zero(poly, ctx);
    }
    else
    {
        ca_poly_fit_length(poly, 1, ctx);
        ca_set(poly->coeffs, x, ctx);
        _ca_poly_set_length(poly, 1, ctx);
    }
}
