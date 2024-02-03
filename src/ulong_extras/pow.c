/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong _n_pow_check(ulong n, ulong exp)
{
    ulong ix, rl, rh;

    if (n == 0)
        FLINT_UNREACHABLE;

    rl = UWORD(1);
    for (ix = 0; ix < exp; ix++)
    {
        umul_ppmm(rh, rl, n, rl);
        if (rh)
            return 0;
    }

    return rl;
}

ulong n_pow(ulong n, ulong exp)
{
    ulong ix;
    ulong res;

    res = UWORD(1);
    for (ix = 0; ix < exp; ix++)
        res *= n;

    return res;
}
