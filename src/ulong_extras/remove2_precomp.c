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
n_remove2_precomp(ulong * n, ulong p, double ppre)
{
    int exp = 0;
    ulong quot, rem = UWORD(0);

    if (p == 2)
    {
        exp = flint_ctz(*n);
        if (exp)
            (*n) >>= exp;

        return exp;
    }

    do
    {
        if ((*n) < p)
            break;
        rem = n_divrem2_precomp(&quot, *n, p, ppre);
        if (rem)
            break;
        exp++;
        (*n) = quot;
    } while (1);

    return exp;
}
