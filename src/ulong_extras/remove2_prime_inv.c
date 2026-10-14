/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

int
n_remove2_prime_inv(ulong * n, ulong p, ulong inv1, ulong inv2)
{
    int val = 0;

    if (p == 2)
    {
        val = flint_ctz(*n);

        if (val)
            (*n) >>= val;
    }

    else
    {
        while (n_divisible_odd_gm(*n, inv1, inv2))
        {
            (*n) *= inv1;
            val ++;
        }
    }

    return val;
}
