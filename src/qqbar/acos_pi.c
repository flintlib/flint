/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "qqbar.h"

int
qqbar_acos_pi(slong * p, ulong * q, const qqbar_t x)
{
    if (qqbar_asin_pi(p, q, x))
    {
        slong a, b, g;
        a = *p;
        b = *q;

        /* 1/2 - a/b */
        a = b - 2 * a;
        b = 2 * b;

        g = n_gcd(FLINT_ABS(a), b);
        if (g != 1)
        {
            a /= g;
            b /= g;
        }

        *p = a;
        *q = b;

        return 1;
    }

    return 0;
}
