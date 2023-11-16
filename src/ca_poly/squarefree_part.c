/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

int
ca_poly_squarefree_part(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
{
    ca_poly_t t;
    int success;

    if (poly->length <= 1)
    {
        ca_poly_one(res, ctx);
        return 1;
    }

    if (poly->length == 2)
        return ca_poly_make_monic(res, poly, ctx);

    ca_poly_init(t, ctx);
    ca_poly_derivative(t, poly, ctx);
    success = ca_poly_gcd(t, poly, t, ctx);
    if (success)
    {
        if (t->length == 1)  /* gcd = 1 */
        {
            success = ca_poly_make_monic(res, poly, ctx);
        }
        else
        {
            success = ca_poly_div(res, poly, t, ctx);
            if (success)
                success = ca_poly_make_monic(res, res, ctx);
        }
    }
    ca_poly_clear(t, ctx);

    return success;
}
