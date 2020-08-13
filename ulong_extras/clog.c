/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_clog(mp_limb_t n, mp_limb_t b)
{
    mp_limb_t r, p, t, phi;

    r = 0;
    p = 1;

    while (1)
    {
        umul_ppmm(phi, t, p, b);

        if (t <= n && !phi)
        {
            r++;
            p = t;
        }
        else
            return r + (p != n);
    }
}
