/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_interval.h"

di_t di_fast_mul(di_t x, di_t y)
{
    di_t res;

    if (x.a > 0 && y.a > 0)
    {
        res.a = x.a * y.a;
        res.b = x.b * y.b;
    }
    else if (x.a > 0 && y.b < 0)
    {
        res.a = x.b * y.a;
        res.b = x.a * y.b;
    }
    else if (x.b < 0 && y.a > 0)
    {
        res.a = x.a * y.b;
        res.b = x.b * y.a;
    }
    else if (x.b < 0 && y.b < 0)
    {
        res.a = x.b * y.b;
        res.b = x.a * y.a;
    }
    else
    {
        double a, b, c, d;

        a = x.a * y.a;
        b = x.a * y.b;
        c = x.b * y.a;
        d = x.b * y.b;

        if (a != a || b != b || c != c || d != d)
        {
            res.a = -D_INF;
            res.b = D_INF;
        }
        else
        {
            res.a = FLINT_MIN(FLINT_MIN(a, b), FLINT_MIN(c, d));
            res.b = FLINT_MAX(FLINT_MAX(a, b), FLINT_MAX(c, d));
        }
    }

    res.a = _di_below(res.a);
    res.b = _di_above(res.b);

    return res;
}

