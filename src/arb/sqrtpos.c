/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_sqrtpos(arb_t z, const arb_t x, slong prec)
{
    if (!arb_is_finite(x))
    {
        if (mag_is_zero(arb_radref(x)) && arf_is_pos_inf(arb_midref(x)))
            arb_pos_inf(z);
        else
            arb_zero_pm_inf(z);
    }
    else if (arb_contains_nonpositive(x))
    {
        arf_t t;

        arf_init(t);

        arf_set_mag(t, arb_radref(x));
        arf_add(t, arb_midref(x), t, MAG_BITS, ARF_RND_CEIL);

        if (arf_sgn(t) <= 0)
        {
            arb_zero(z);
        }
        else
        {
            arf_sqrt(t, t, MAG_BITS, ARF_RND_CEIL);
            arf_mul_2exp_si(t, t, -1);
            arf_set(arb_midref(z), t);
            arf_get_mag(arb_radref(z), t);
        }

        arf_clear(t);
    }
    else
    {
        arb_sqrt(z, x, prec);
    }

    arb_nonnegative_part(z, z);
}

