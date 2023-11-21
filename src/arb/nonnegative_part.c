/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_nonnegative_part(arb_t res, const arb_t x)
{
    if (!arb_contains_negative(x))
    {
        arb_set(res, x);
    }
    else if (!arb_is_finite(x))
    {
        arb_indeterminate(res);
    }
    else
    {
        arf_t t;
        arf_init(t);

        arf_set_mag(t, arb_radref(x));
        arf_add(arb_midref(res), arb_midref(x), t, MAG_BITS, ARF_RND_CEIL);

        if (arf_sgn(arb_midref(res)) <= 0)
        {
            arf_zero(arb_midref(res));
            mag_zero(arb_radref(res));
        }
        else
        {
            arf_mul_2exp_si(arb_midref(res), arb_midref(res), -1);
            arf_get_mag(arb_radref(res), arb_midref(res));
            /* needed since arf_get_mag is inexact */
            arf_set_mag(arb_midref(res), arb_radref(res));
        }

        arf_clear(t);
    }
}

