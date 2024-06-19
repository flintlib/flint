/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong
n_powmod_ui_precomp(ulong a, ulong exp, ulong n, double npre)
{
    ulong x, y;

    if (n == UWORD(1))
        return WORD(0);

    x = UWORD(1);
    y = a;

    while (exp)
    {
        if (exp & WORD(1))
            x = n_mulmod_precomp(x, y, n, npre);
        exp >>= 1;
        if (exp)
            y = n_mulmod_precomp(y, y, n, npre);
    }

    return x;
}

ulong
n_powmod_precomp(ulong a, slong exp, ulong n, double npre)
{
    if (exp < 0)
    {
        a = n_invmod(a, n);
        exp = -exp;
    }

    return n_powmod_ui_precomp(a, exp, n, npre);
}
