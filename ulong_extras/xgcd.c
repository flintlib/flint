/*
    Copyright (C) 2009, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

ulong
n_xgcd(ulong * a, ulong * b, ulong x, ulong y)
{
    slong u1, u2, v1, v2, t1, t2;
    ulong u3, v3, quot, rem, d;

    FLINT_ASSERT(x >= y);

    u1 = v2 = 1;
    u2 = v1 = 0;
    u3 = x;
    v3 = y;

    /* x and y both have top bit set */
    if ((slong) (x & y) < WORD(0))
    {
        d = u3 - v3;
        t2 = v2;
        t1 = u2;
        u2 = u1 - u2;
        u1 = t1;
        u3 = v3;
        v2 = v1 - v2;
        v1 = t2;
        v3 = d;
    }

    /* second value has second msb set */
    while ((slong) (v3 << 1) < WORD(0))
    {
        d = u3 - v3;
        if (d < v3)             /* quot = 1 */
        {
            t2 = v2;
            t1 = u2;
            u2 = u1 - u2;
            u1 = t1;
            u3 = v3;
            v2 = v1 - v2;
            v1 = t2;
            v3 = d;
        }
        else if (d < (v3 << 1)) /* quot = 2 */
        {
            t1 = u2;
            u2 = u1 - (u2 << 1);
            u1 = t1;
            u3 = v3;
            t2 = v2;
            v2 = v1 - (v2 << 1);
            v1 = t2;
            v3 = d - u3;
        }
        else                    /* quot = 3 */
        {
            t1 = u2;
            u2 = u1 - 3 * u2;
            u1 = t1;
            u3 = v3;
            t2 = v2;
            v2 = v1 - 3 * v2;
            v1 = t2;
            v3 = d - (u3 << 1);
        }
    }

    while (v3)
    {
        d = u3 - v3;

        /* overflow not possible, top 2 bits of v3 not set */
        if (u3 < (v3 << 2))
        {
            if (d < v3)         /* quot = 1 */
            {
                t2 = v2;
                t1 = u2;
                u2 = u1 - u2;
                u1 = t1;
                u3 = v3;
                v2 = v1 - v2;
                v1 = t2;
                v3 = d;
            }
            else if (d < (v3 << 1)) /* quot = 2 */
            {
                t1 = u2;
                u2 = u1 - (u2 << 1);
                u1 = t1;
                u3 = v3;
                t2 = v2;
                v2 = v1 - (v2 << 1);
                v1 = t2;
                v3 = d - u3;
            }
            else                /* quot = 3 */
            {
                t1 = u2;
                u2 = u1 - 3 * u2;
                u1 = t1;
                u3 = v3;
                t2 = v2;
                v2 = v1 - 3 * v2;
                v1 = t2;
                v3 = d - (u3 << 1);
            }
        }
        else
        {
            quot = u3 / v3;
            rem = u3 - v3 * quot;
            t1 = u2;
            u2 = u1 - quot * u2;
            u1 = t1;
            u3 = v3;
            t2 = v2;
            v2 = v1 - quot * v2;
            v1 = t2;
            v3 = rem;
        }
    }

    /* Remarkably, |u1| < x/2, thus comparison with 0 is valid */
    if (u1 <= WORD(0))
    {
        u1 += y;
        v1 -= x;
    }

    *a = u1;
    *b = -v1;

    return u3;
}
