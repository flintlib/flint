/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

int
n_remove(ulong * n, ulong p)
{
    int exp, i;
    ulong powp[6];
    ulong quot, rem;

    if (p == 2)
    {
        exp = flint_ctz(*n);
        if (exp)
            (*n) >>= exp;

        return exp;
    }

    powp[0] = p;

    for (i = 0; ; i++)
    {
        if ((*n) < powp[i])
            break;
        quot = (*n) / powp[i];
        rem = (*n) - quot * powp[i];
        if (rem != UWORD(0))
            break;
        powp[i + 1] = powp[i] * powp[i];
        (*n) = quot;
    }

    exp = (1 << i) - 1;

    while (i > 0)
    {
        i--;
        if ((*n) < powp[i])
            continue;
        quot = (*n) / powp[i];
        rem = (*n) - quot * powp[i];
        if (rem == UWORD(0))
        {
            exp += (UWORD(1) << i);
            (*n) = quot;
        }
    }

    return exp;
}
