/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_contains_int(const arb_t x)
{
    if (arf_is_int(arb_midref(x)))
    {
        return 1;
    }
    else if (!arb_is_finite(x))
    {
        return arb_contains_zero(x);
    }
    else if (arb_is_exact(x))
    {
        return 0;
    }
    else if (mag_cmp_2exp_si(arb_radref(x), -1) >= 0)  /* radius >= 1/2 */
    {
        return 1;
    }
    else
    {
        /* radius is < 1/2, so it's enough to test the two integers
            bracketing the midpoint */
        arf_t t;
        int res;
        arf_init(t);
        arf_floor(t, arb_midref(x));

        res = arb_contains_arf(x, t);
        if (!res)
        {
            arf_ceil(t, arb_midref(x));
            res = arb_contains_arf(x, t);
        }

        arf_clear(t);
        return res;
    }
}

