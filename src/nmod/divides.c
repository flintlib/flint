/*
    Copyright (C) 2019 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"

int nmod_divides(ulong * a, ulong b, ulong c, nmod_t mod)
{
    int success;
    ulong g, x, y, q;

    if (c == 0)
    {
        if (b == 0)
        {
            *a = 0;
            return 1;
        }
        else
        {
            *a = 0;
            return 0;
        }
    }
    else if (b == 0)
    {
        *a = 0;
        return 1;
    }

    /* solve g = c*(-x) + n*y where g = gcd(c, n) */
    g = n_xgcd(&y, &x, mod.n, c);
    success = (b % g == 0);
    if (success)
    {
        q = b / g;
        *a = nmod_mul(q, nmod_neg(x, mod), mod);
    }

    return success;
}
