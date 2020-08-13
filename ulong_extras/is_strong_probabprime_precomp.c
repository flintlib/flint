/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Dana Jacobsen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_is_strong_probabprime_precomp(mp_limb_t n, double npre, mp_limb_t a,
                                mp_limb_t d)
{
    mp_limb_t t = d;
    mp_limb_t y;

    /* Map large base to range 2 ... n - 1 */
    if (a >= n)  
       a = n_mod2_precomp(a, n, npre);

    if ((a <= 1) || (a == n-1))  return 1;

    y = n_powmod_ui_precomp(a, t, n, npre);

    if (y == UWORD(1))
        return 1;
    t <<= 1;

    while ((t != n - 1) && (y != n - 1))
    {
        y = n_mulmod_precomp(y, y, n, npre);
        t <<= 1;
    }

    return (y == n - 1);
}
