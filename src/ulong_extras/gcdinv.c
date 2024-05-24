/*
    Copyright (C) 2009, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

#if N_GCDEXT_METHOD == 0
nn_pair_t
_n_gcdinv(ulong x, ulong y)
{
    nn_pair_t ret;
    slong v1, v2, t2;
    ulong d, r, quot, rem;

    FLINT_ASSERT(y > x);

    v1 = 0;
    v2 = 1;
    r = x;
    x = y;

    /* y and x both have top bit set */
    if ((slong) (x & r) < 0)
    {
        d = x - r;
        t2 = v2;
        x = r;
        v2 = v1 - v2;
        v1 = t2;
        r = d;
    }

    /* second value has second msb set */
    while ((slong) (r << 1) < 0)
    {
        d = x - r;
        if (d < r)              /* quot = 1 */
        {
            t2 = v2;
            x = r;
            v2 = v1 - v2;
            v1 = t2;
            r = d;
        }
        else if (d < (r << 1))  /* quot = 2 */
        {
            x = r;
            t2 = v2;
            v2 = v1 - (v2 << 1);
            v1 = t2;
            r = d - x;
        }
        else                    /* quot = 3 */
        {
            x = r;
            t2 = v2;
            v2 = v1 - 3 * v2;
            v1 = t2;
            r = d - (x << 1);
        }
    }

    while (r)
    {
        /* overflow not possible due to top 2 bits of r not being set */
        if (x < (r << 2))       /* if quot < 4 */
        {
            d = x - r;
            if (d < r)          /* quot = 1 */
            {
                t2 = v2;
                x = r;
                v2 = v1 - v2;
                v1 = t2;
                r = d;
            }
            else if (d < (r << 1))  /* quot = 2 */
            {
                x = r;
                t2 = v2;
                v2 = v1 - (v2 << 1);
                v1 = t2;
                r = d - x;
            }
            else                /* quot = 3 */
            {
                x = r;
                t2 = v2;
                v2 = v1 - 3 * v2;
                v1 = t2;
                r = d - (x << 1);
            }
        }
        else
        {
            quot = x / r;
            rem = x - r * quot;
            x = r;
            t2 = v2;
            v2 = v1 - quot * v2;
            v1 = t2;
            r = rem;
        }
    }

    if (v1 < WORD(0))
        v1 += y;

    ret.m0 = v1;
    ret.m1 = x;

    return ret;
}
#elif N_GCDEXT_METHOD == 1
nn_pair_t
_n_gcdinv(ulong x, ulong y)
{
    nn_pair_t ret;
    slong v1, v2, t2;
    ulong d, r, quot, rem;

    FLINT_ASSERT(y > x);

    v1 = 0;
    v2 = 1;
    r = x;
    x = y;

    while (r)
    {
        quot = x / r;
        rem = x - r * quot;
        x = r;
        t2 = v2;
        v2 = v1 - quot * v2;
        v1 = t2;
        r = rem;
    }

    if (v1 < WORD(0))
        v1 += y;

    ret.m0 = v1;
    ret.m1 = x;

    return ret;
}
#else
# error
#endif
