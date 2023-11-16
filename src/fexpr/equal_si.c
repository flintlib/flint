/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

int
fexpr_equal_si(const fexpr_t expr, slong c)
{
    if (c >= FEXPR_COEFF_MIN && c <= FEXPR_COEFF_MAX)
        return expr->data[0] == (c << FEXPR_TYPE_BITS);
    else if (c > 0)
        return (expr->data[0] == (FEXPR_TYPE_BIG_INT_POS | (2 << FEXPR_TYPE_BITS)) && expr->data[1] == c);
    else
        return (expr->data[0] == (FEXPR_TYPE_BIG_INT_NEG | (2 << FEXPR_TYPE_BITS)) && expr->data[1] == (-(ulong) c));
}
