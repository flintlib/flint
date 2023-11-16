/*
    Copyright (C) 2022 Erik Postma

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_nonnegative_abs(arb_t y, const arb_t x)
{
    if(arb_is_finite(x) && arb_contains_zero(x))
    {
        /* We need to round to MAG_BITS down below, anyway. */
        arb_get_abs_ubound_arf(arb_midref(y), x, MAG_BITS+1);

        /* Now t := arb_midref(y) is the upper bound of the interval x; we
           need to set arb_midref(y) and arb_radref(y) to (approximations
           of) t/2.  */
        arf_mul_2exp_si(arb_midref(y), arb_midref(y), -1);
        arf_get_mag(arb_radref(y), arb_midref(y));
        /* The above is inexact (rounding up), so we need to update
           arb_midref(y) to match arb_radref(y) again. (That is exact.) */
        arf_set_mag(arb_midref(y), arb_radref(y));
    }
    else
    {
        arf_abs(arb_midref(y), arb_midref(x));
        mag_set(arb_radref(y), arb_radref(x));
    }
}
