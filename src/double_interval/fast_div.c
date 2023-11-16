/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_interval.h"

di_t di_fast_div(di_t x, di_t y)
{
    di_t res;

    if (y.a > 0)
    {
        if (x.a >= 0)
        {
            res.a = x.a / y.b;
            res.b = x.b / y.a;
        }
        else if (x.b <= 0)
        {
            res.a = x.a / y.a;
            res.b = x.b / y.b;
        }
        else
        {
            res.a = x.a / y.a;
            res.b = x.b / y.a;
        }
    }
    else if (y.b < 0)
    {
        if (x.a >= 0)
        {
            res.a = x.b / y.b;
            res.b = x.a / y.a;
        }
        else if (x.b <= 0)
        {
            res.a = x.b / y.a;
            res.b = x.a / y.b;
        }
        else
        {
            res.a = x.b / y.b;
            res.b = x.a / y.b;
        }
    }
    else
    {
        res.a = -D_INF;
        res.b = D_INF;
    }

    res.a = _di_below(res.a);
    res.b = _di_above(res.b);

    return res;
}
