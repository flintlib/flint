/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

/* todo: bounds checking? */
void
fexpr_arg(fexpr_t res, const fexpr_t expr, slong i)
{
    const ulong * data;
    slong j, size;
    ulong type = FEXPR_TYPE(expr->data[0]);

    if (FEXPR_TYPE_CALL0 <= type && type <= FEXPR_TYPE_CALL4)
    {
        data = expr->data + FEXPR_HEADER_SIZE;
        data += FEXPR_SIZE(data[0]);  /* skip f */

        for (j = 0; j < i; j++)
            data += FEXPR_SIZE(data[0]);   /* jump ahead */

        size = FEXPR_SIZE(data[0]);
        fexpr_fit_size(res, size);
        flint_mpn_copyi(res->data, data, size);
    }
    else if (type == FEXPR_TYPE_CALLN)
    {
        data = expr->data + expr->data[3 + i / 4];

        for (j = 0; j < i % 4; j++)
            data += FEXPR_SIZE(data[0]);

        size = FEXPR_SIZE(data[0]);
        fexpr_fit_size(res, size);
        flint_mpn_copyi(res->data, data, size);
    }
    else
    {
        flint_throw(FLINT_ERROR, "fexpr_arg: a non-atomic expression is required\n");
    }
}

void
fexpr_view_arg(fexpr_t res, const fexpr_t expr, slong i)
{
    const ulong * data;
    slong j;
    ulong type = FEXPR_TYPE(expr->data[0]);

    if (FEXPR_TYPE_CALL0 <= type && type <= FEXPR_TYPE_CALL4)
    {
        data = expr->data + FEXPR_HEADER_SIZE;
        data += FEXPR_SIZE(data[0]);  /* skip f */

        for (j = 0; j < i; j++)
            data += FEXPR_SIZE(data[0]);   /* jump ahead */

        res->data = (ulong *) data;
        res->alloc = 0;
    }
    else if (type == FEXPR_TYPE_CALLN)
    {
        data = expr->data + expr->data[3 + i / 4];

        for (j = 0; j < i % 4; j++)
            data += FEXPR_SIZE(data[0]);

        res->data = (ulong *) data;
        res->alloc = 0;
    }
    else
    {
        flint_throw(FLINT_ERROR, "fexpr_view_arg: a non-atomic expression is required\n");
    }
}
