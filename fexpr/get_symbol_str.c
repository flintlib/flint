/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

char * fexpr_get_symbol_str(const fexpr_t expr)
{
    char * res;
    slong i, len;
    ulong head = expr->data[0];

    if (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_SYMBOL)
    {
        res = flint_malloc(FEXPR_SMALL_SYMBOL_LEN + 1);
        res[FEXPR_SMALL_SYMBOL_LEN] = '\0';

        for (i = 0; i < FEXPR_SMALL_SYMBOL_LEN; i++)
        {
            res[i] = (head >> ((i + 1) * 8));
            if (res[i] == '\0')
                break;
        }
    }
    else if (FEXPR_TYPE(head) == FEXPR_TYPE_BIG_SYMBOL)
    {
        len = strlen((const char *) (expr->data + 1));
        res = flint_malloc(len + 1);
        memcpy(res, (const char *) (expr->data + 1), len + 1);
    }
    else
    {
        flint_printf("fexpr_get_symbol_str: a symbol is required\n");
        flint_abort();
    }

    return res;
}
