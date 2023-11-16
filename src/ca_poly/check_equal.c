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
_ca_poly_check_equal(ca_srcptr poly1, slong len1,
        ca_srcptr poly2, slong len2, ca_ctx_t ctx)
{
    truth_t eq, res;
    slong i;

    res = T_TRUE;

    for (i = 0; i < len2; i++)
    {
        eq = ca_check_equal(poly1 + i, poly2 + i, ctx);

        if (eq == T_FALSE)
            return T_FALSE;
        if (eq == T_UNKNOWN)
            res = T_UNKNOWN;
    }

    for (i = len2; i < len1; i++)
    {
        eq = ca_check_is_zero(poly1 + i, ctx);

        if (eq == T_FALSE)
            return T_FALSE;
        if (eq == T_UNKNOWN)
            res = T_UNKNOWN;
    }

    return res;
}

truth_t
ca_poly_check_equal(const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;

    if (len1 >= len2)
        return _ca_poly_check_equal(poly1->coeffs, len1, poly2->coeffs, len2, ctx);
    else
        return _ca_poly_check_equal(poly2->coeffs, len2, poly1->coeffs, len1, ctx);
}
