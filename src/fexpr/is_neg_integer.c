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
fexpr_is_neg_integer(const fexpr_t expr)
{
    ulong head = expr->data[0];

    if (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_INT)
        return (((slong) head) >> FEXPR_TYPE_BITS) < 0;

    if (FEXPR_TYPE(head) == FEXPR_TYPE_BIG_INT_NEG)
        return 1;

    return 0;
}
