/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"

ulong n_mulmod_precomp(ulong a, ulong b, ulong n, double npre)
{
    ulong quot;
    slong rem;

    quot = (ulong) ((double) a * (double) b * npre);
    rem  = a * b - quot * n;
    if (rem < 0) 
    {
        rem += n;
        if (rem < 0) return rem + n;
    }
    else if (rem >= n) return rem - n;
    return rem;
}
