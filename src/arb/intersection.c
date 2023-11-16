/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
arb_intersection(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    arf_t left, right, t, xr, yr;
    int result;

    if (arf_is_nan(arb_midref(x)) || arf_is_nan(arb_midref(y)))
    {
        arb_indeterminate(z);
        return 1;
    }

    if (mag_is_inf(arb_radref(x)) && mag_is_inf(arb_radref(y)))
    {
        arb_zero_pm_inf(z);
        return 1;
    }

    result = arb_overlaps(x, y);

    if (result)
    {
        arf_init(left);
        arf_init(right);
        arf_init(t);

        arf_init_set_mag_shallow(xr, arb_radref(x));
        arf_init_set_mag_shallow(yr, arb_radref(y));

        arf_sub(left, arb_midref(x), xr, prec, ARF_RND_FLOOR);
        arf_sub(t, arb_midref(y), yr, prec, ARF_RND_FLOOR);
        arf_max(left, left, t);

        arf_add(right, arb_midref(x), xr, prec, ARF_RND_CEIL);
        arf_add(t, arb_midref(y), yr, prec, ARF_RND_CEIL);
        arf_min(right, right, t);

        arb_set_interval_arf(z, left, right, prec);

        arf_clear(left);
        arf_clear(right);
        arf_clear(t);
    }

    return result;
}
