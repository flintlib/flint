/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

ulong n_flog(ulong n, ulong b)
{
    ulong r, p, phi;

    r = 0;
    p = 1;

    while (1)
    {
        umul_ppmm(phi, p, p, b);

        if (p <= n && !phi)
            r++;
        else
            return r;
    }
}
