/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int
ca_qqbar_acos_pi(slong * p, ulong * q, const ca_qqbar_t x)
{
    if (ca_qqbar_asin_pi(p, q, x))
    {
        slong a, b;
        a = *p;
        b = *q;

        /* 1/2 - a/b */
        if (b % 2 == 0)
        {
            *p = b / 2 - a;
            *q = b;
        }
        else
        {
            *p = b - 2 * a;
            *q = 2 * b;
        }

        return 1;
    }

    return 0;
}
