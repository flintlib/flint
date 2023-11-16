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
ca_poly_is_proper(const ca_poly_t poly, ca_ctx_t ctx)
{
    slong i, len;

    len = poly->length;

    for (i = 0; i < len; i++)
        if (CA_IS_SPECIAL(poly->coeffs + i))
            return 0;

    if (len >= 1)
        if (ca_check_is_zero(poly->coeffs + len - 1, ctx) != T_FALSE)
            return 0;

    return 1;
}
