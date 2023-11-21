/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

int
mag_cmp(const mag_t x, const mag_t y)
{
    int c;

    if (mag_equal(x, y))
        return 0;

    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_zero(x)) return -1;
        if (mag_is_zero(y)) return 1;
        if (mag_is_inf(x)) return 1;
        if (mag_is_inf(y)) return -1;
    }

    c = fmpz_cmp(MAG_EXPREF(x), MAG_EXPREF(y));

    if (c == 0)
        return (MAG_MAN(x) < MAG_MAN(y)) ? -1 : 1;

    return (c < 0) ? -1 : 1;
}

