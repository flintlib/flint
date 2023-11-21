/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_get_interval_arf(arf_t a, arf_t b, const arb_t x, slong prec)
{
    arf_t r;
    arf_init_set_mag_shallow(r, arb_radref(x));
    arf_sub(a, arb_midref(x), r, prec, ARF_RND_FLOOR);
    arf_add(b, arb_midref(x), r, prec, ARF_RND_CEIL);
}

