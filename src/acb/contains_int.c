/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

int
acb_contains_int(const acb_t x)
{
    if (!arb_contains_zero(acb_imagref(x)))
        return 0;

    return arb_contains_int(acb_realref(x));
}
