/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_hypot(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_zero(y))
    {
        mag_set(z, x);
    }
    else if (mag_is_zero(x))
    {
        mag_set(z, y);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_mul(t, x, x);
        mag_addmul(t, y, y);
        mag_sqrt(z, t);
        mag_clear(t);
    }
}

