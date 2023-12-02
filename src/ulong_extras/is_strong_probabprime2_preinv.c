/*
    Copyright (C) 2008, Peter Shrimpton
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Dana Jacobsen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

int
n_is_strong_probabprime2_preinv(mp_limb_t n, mp_limb_t ninv, mp_limb_t a,
                                mp_limb_t d)
{
    mp_limb_t t = d;
    mp_limb_t y;

    FLINT_ASSERT(a < n);

    if (a <= 1 || a == n - 1)
        return 1;

    y = n_powmod2_ui_preinv(a, t, n, ninv);

    if (y == UWORD(1))
        return 1;
    t <<= 1;

    while ((t != n - 1) && (y != n - 1))
    {
        y = n_mulmod2_preinv(y, y, n, ninv);
        t <<= 1;
    }

    return (y == n - 1);
}
