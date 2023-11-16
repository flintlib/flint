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
arb_get_abs_ubound_arf(arf_t u, const arb_t x, slong prec)
{
    arf_t t;
    arf_init_set_mag_shallow(t, arb_radref(x));

    if (arf_sgn(arb_midref(x)) < 0)
        arf_sub(u, arb_midref(x), t, prec, ARF_RND_UP);
    else
        arf_add(u, arb_midref(x), t, prec, ARF_RND_UP);

    arf_abs(u, u);
}
