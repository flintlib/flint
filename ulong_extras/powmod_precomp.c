/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t
n_powmod_ui_precomp(mp_limb_t a, mp_limb_t exp, mp_limb_t n, double npre)
{
    mp_limb_t x, y;

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

mp_limb_t
n_powmod_precomp(mp_limb_t a, mp_limb_signed_t exp, mp_limb_t n, double npre)
{
    if (exp < 0)
    {
        a = n_invmod(a, n);
        exp = -exp;
    }

    return n_powmod_ui_precomp(a, exp, n, npre);
}
