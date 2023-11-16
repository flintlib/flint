/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

ulong
fexpr_hash(const fexpr_t expr)
{
    ulong head, hash;
    slong i, size;

    hash = head = expr->data[0];
    size = FEXPR_SIZE(head);

    for (i = 1; i < size; i++)
        hash = expr->data[i] * 1000003 + hash;

    return hash;
}
