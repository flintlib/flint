/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_ceil(arb_t z, const arb_t x, slong prec)
{
    if (!arb_is_finite(x))
    {
        arb_indeterminate(z);
    }
    else if (arb_is_exact(x))
    {
        arf_ceil(arb_midref(z), arb_midref(x));
        mag_zero(arb_radref(z));
        arb_set_round(z, z, prec);
    }
    else
    {
        arf_t a, b;
        arf_init(a);
        arf_init(b);

        arb_get_interval_arf(a, b, x, prec);
        arf_ceil(a, a);
        arf_ceil(b, b);
        arb_set_interval_arf(z, a, b, prec);

        arf_clear(a);
        arf_clear(b);
    }
}

