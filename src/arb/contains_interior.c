/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_contains_interior(const arb_t x, const arb_t y)
{
    arf_t t;
    arf_t u;
    arf_t xr, yr;
    arf_struct tmp[4];
    int left_ok, right_ok;

    if (arf_is_nan(arb_midref(x)) || arb_is_exact(x) || !arb_is_finite(y))
        return 0;

    arf_init(t);
    arf_init(u);

    arf_init_set_mag_shallow(xr, arb_radref(x));
    arf_init_set_mag_shallow(yr, arb_radref(y));

    /* fast check */
    arf_sub(t, arb_midref(x), xr, MAG_BITS, ARF_RND_CEIL);
    arf_sub(u, arb_midref(y), yr, MAG_BITS, ARF_RND_FLOOR);
    left_ok = arf_cmp(t, u) < 0;

    /* exact check */
    if (!left_ok)
    {
        arf_init_set_shallow(tmp + 0, arb_midref(x));
        arf_init_neg_mag_shallow(tmp + 1, arb_radref(x));
        arf_init_neg_shallow(tmp + 2, arb_midref(y));
        arf_init_set_mag_shallow(tmp + 3, arb_radref(y));

        arf_sum(t, tmp, 4, MAG_BITS, ARF_RND_DOWN);
        left_ok = arf_sgn(t) < 0;
    }

    /* fast check */
    arf_add(t, arb_midref(x), xr, MAG_BITS, ARF_RND_FLOOR);
    arf_add(u, arb_midref(y), yr, MAG_BITS, ARF_RND_CEIL);
    right_ok = (arf_cmp(t, u) > 0);

    /* exact check */
    if (!right_ok)
    {
        arf_init_set_shallow(tmp + 0, arb_midref(x));
        arf_init_set_mag_shallow(tmp + 1, arb_radref(x));
        arf_init_neg_shallow(tmp + 2, arb_midref(y));
        arf_init_neg_mag_shallow(tmp + 3, arb_radref(y));

        arf_sum(t, tmp, 4, MAG_BITS, ARF_RND_DOWN);
        right_ok = arf_sgn(t) > 0;
    }

    arf_clear(t);
    arf_clear(u);

    return left_ok && right_ok;
}

