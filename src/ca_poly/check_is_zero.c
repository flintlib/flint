/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

truth_t
ca_poly_check_is_zero(const ca_poly_t poly, ca_ctx_t ctx)
{
    ca_t t;
    truth_t res;

    if (poly->length == 0)
        return T_TRUE;

    ca_init(t, ctx);
    res = _ca_poly_check_equal(poly->coeffs, poly->length, t, 1, ctx);
    ca_clear(t, ctx);
    return res;
}
