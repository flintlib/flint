/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_interval.h"

di_t di_fast_sqr(di_t x)
{
    di_t res;

    if (x.a >= 0)
    {
        res.a = x.a * x.a;
        res.b = x.b * x.b;
    }
    else if (x.b <= 0)
    {
        res.a = x.b * x.b;
        res.b = x.a * x.a;
    }
    else
    {
        res.a = 0.0;
        res.b = FLINT_MAX(x.a * x.a, x.b * x.b);
    }

    if (res.a != 0.0)
        res.a = _di_below(res.a);

    res.b = _di_above(res.b);

    return res;
}

