/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"

ulong n_sqrtmod(ulong a, ulong p)
{
    nmod_t mod;
    ulong r;

    if (a <= 1)
        return a;

    /* skip nmod_init, which costs as much as the search */
    if (p < NMOD_SQRT_SMALL_CUTOFF && (p & 1))
    {
        _nmod_sqrt_small(&r, a, p);
        return r;
    }

    nmod_init(&mod, p);
    nmod_sqrt(&r, a, mod);
    return r;
}
