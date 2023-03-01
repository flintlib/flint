/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int
acb_get_unique_fmpz(fmpz_t z, const acb_t x)
{
    if (!arb_contains_zero(acb_imagref(x)))
        return 0;

    return arb_get_unique_fmpz(z, acb_realref(x));
}

