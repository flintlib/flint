/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_sgn(arb_t res, const arb_t x)
{
    if (arb_is_zero(x))
    {
        arb_zero(res);
    }
    else if (arb_contains_zero(x))
    {
        arf_zero(arb_midref(res));
        mag_one(arb_radref(res));
    }
    else
    {
        arb_set_si(res, arf_sgn(arb_midref(x)));
    }
}

