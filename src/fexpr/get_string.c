/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"

#ifdef __GNUC__
# define memcpy __builtin_memcpy
# define strlen __builtin_strlen
#else
# include <string.h>
#endif

char * fexpr_get_string(const fexpr_t expr)
{
    char * res;
    slong i, len;
    ulong head = expr->data[0];

    if (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_STRING)
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
    else if (FEXPR_TYPE(head) == FEXPR_TYPE_BIG_STRING)
    {
        len = strlen((const char *) (expr->data + 1));
        res = flint_malloc(len + 1);
        memcpy(res, (const char *) (expr->data + 1), len + 1);
    }
    else
    {
        flint_throw(FLINT_ERROR, "fexpr_get_string: a string is required\n");
    }

    return res;
}
