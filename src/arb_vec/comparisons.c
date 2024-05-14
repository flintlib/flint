/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_vec.h"

int _arb_vec_is_zero(arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!arb_is_zero(vec + i))
            return 0;
    return 1;
}

int _arb_vec_is_finite(arb_srcptr x, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        if (!arb_is_finite(x + i))
            return 0;

    return 1;
}

int _arb_vec_equal(arb_srcptr vec1, arb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!arb_equal(vec1 + i, vec2 + i))
            return 0;
    }
    return 1;
}

int _arb_vec_overlaps(arb_srcptr vec1, arb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!arb_overlaps(vec1 + i, vec2 + i))
            return 0;
    }

    return 1;
}

int _arb_vec_contains(arb_srcptr vec1, arb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!arb_contains(vec1 + i, vec2 + i))
            return 0;
    }

    return 1;
}
