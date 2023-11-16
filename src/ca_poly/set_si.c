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
ca_poly_set_si(ca_poly_t poly, slong x, ca_ctx_t ctx)
{
    if (x == 0)
    {
        ca_poly_zero(poly, ctx);
    }
    else
    {
        ca_poly_fit_length(poly, 1, ctx);
        ca_set_si(poly->coeffs, x, ctx);
        _ca_poly_set_length(poly, 1, ctx);
    }
}
