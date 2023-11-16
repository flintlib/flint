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
fexpr_cmp_fast(const fexpr_t a, const fexpr_t b)
{

    ulong ha, hb;
    slong sa, sb;
    slong i;

    ha = a->data[0];
    hb = b->data[0];

    if (ha != hb)
        return (ha > hb) ? 1 : -1;

    sa = FEXPR_SIZE(ha);
    sb = FEXPR_SIZE(hb);

    if (sa != sb)
        return 0;

    for (i = 1; i < sa; i++)
    {
        ha = a->data[i];
        hb = b->data[i];
        if (ha != hb)
            return (ha > hb) ? 1 : -1;
    }

    return 0;
}
