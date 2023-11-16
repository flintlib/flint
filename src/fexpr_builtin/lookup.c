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
# define strcmp __builtin_strcmp
#else
# include <string.h>
#endif

slong
fexpr_builtin_lookup(const char * s)
{
    slong a, mid, b;
    int cmp;

    a = 0;
    b = FEXPR_BUILTIN_LENGTH - 1;

    while (a <= b)
    {
        mid = (a + b) / 2;
        cmp = strcmp(fexpr_builtin_table[mid].string, s);

        if (cmp == 0)
            return mid;
        else if (cmp > 0)
            b = mid - 1;
        else
            a = mid + 1;
    }

    return -1;
}
