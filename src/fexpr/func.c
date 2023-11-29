/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

void
fexpr_func(fexpr_t res, const fexpr_t expr)
{
    ulong type = FEXPR_TYPE(expr->data[0]);
    slong size;
    const ulong * data;

    if (FEXPR_TYPE_CALL0 <= type && type <= FEXPR_TYPE_CALL4)
    {
        data = expr->data + FEXPR_HEADER_SIZE;
    }
    else if (type == FEXPR_TYPE_CALLN)
    {
        data = expr->data + expr->data[2];
    }
    else
    {
        flint_throw(FLINT_ERROR, "fexpr_func: a non-atomic expression is required\n");
    }

    size = FEXPR_SIZE(data[0]);
    fexpr_fit_size(res, size);
    flint_mpn_copyi(res->data, data, size);
}

void
fexpr_view_func(fexpr_t res, const fexpr_t expr)
{
    ulong type = FEXPR_TYPE(expr->data[0]);
    const ulong * data;

    if (FEXPR_TYPE_CALL0 <= type && type <= FEXPR_TYPE_CALL4)
    {
        data = expr->data + FEXPR_HEADER_SIZE;
    }
    else if (type == FEXPR_TYPE_CALLN)
    {
        data = expr->data + expr->data[2];
    }
    else
    {
        flint_throw(FLINT_ERROR, "fexpr_view_func: a non-atomic expression is required\n");
    }

    res->data = (ulong *) data;
    res->alloc = 0;
}
