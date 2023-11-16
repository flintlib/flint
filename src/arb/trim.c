/*
    Copyright (C) 2013-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define TRIM_PADDING 16

void
arb_trim(arb_t y, const arb_t x)
{
    if (mag_is_zero(arb_radref(x)) || arf_is_special(arb_midref(x)))
    {
        arb_set(y, x);
    }
    else if (mag_is_special(arb_radref(x)))
    {
        /* midpoint must be finite, so set to 0 +/- inf */
        arb_zero_pm_inf(y);
    }
    else
    {
        slong bits, accuracy;

        bits = arb_bits(x);
        accuracy = arb_rel_accuracy_bits(x);

        if (accuracy < -TRIM_PADDING)
        {
            mag_t t;

            /* set to 0 +/- r */
            mag_init(t);
            arf_get_mag(t, arb_midref(x));
            mag_add(arb_radref(y), t, arb_radref(x));
            mag_clear(t);

            arf_zero(arb_midref(y));
        }
        else if (accuracy < bits - TRIM_PADDING)
        {
            arb_set_round(y, x, FLINT_MAX(0, accuracy) + TRIM_PADDING);
        }
        else
        {
            arb_set(y, x);
        }
    }
}

