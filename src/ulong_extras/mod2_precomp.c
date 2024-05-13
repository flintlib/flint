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

ulong
n_mod2_precomp(ulong a, ulong n, double npre)
{
    ulong quot;
    slong rem;

    if (a < n)
        return a;
    if ((slong) n < WORD(0))
        return a - n;

    if (n == 1)
    {
        quot = a;
        rem = 0;
    } else
    {
        quot = (ulong) ((double) a * npre);
        rem  = a - quot * n;
    }

    if (rem < (slong) (-n))
        quot -= (ulong) ((double) (-rem) * npre);
    else if (rem >= (slong) n)
        quot += (ulong) ((double) rem * npre);
    else if (rem < WORD(0))
        return rem + n;
    else
        return rem;

    rem = a - quot * n;
    if (rem >= (slong) n)
        return rem - n;
    else if (rem < WORD(0))
        return rem + n;
    else
        return rem;
}
