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
ca_poly_make_monic(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
{
    if (poly->length == 0)
    {
        ca_poly_zero(res, ctx);
        return 0;
    }
    else
    {
        if (ca_check_is_one(poly->coeffs + poly->length - 1, ctx) == T_TRUE)
        {
            ca_poly_set(res, poly, ctx);
            ca_one(res->coeffs + res->length - 1, ctx);
            return 1;
        }

        if (ca_check_is_neg_one(poly->coeffs + poly->length - 1, ctx) == T_TRUE)
        {
            ca_poly_neg(res, poly, ctx);
            ca_one(res->coeffs + res->length - 1, ctx);
            return 1;
        }

        ca_poly_set(res, poly, ctx);
        ca_inv(res->coeffs + res->length - 1, res->coeffs + res->length - 1, ctx);
        if (CA_IS_SPECIAL(res->coeffs + res->length - 1))
        {
            return 0;
        }
        else
        {
            _ca_vec_scalar_mul_ca(res->coeffs, res->coeffs, res->length - 1, res->coeffs + res->length - 1, ctx);
            ca_one(res->coeffs + res->length - 1, ctx);
            return 1;
        }
    }

}
