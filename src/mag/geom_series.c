/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_geom_series(mag_t res, const mag_t x, ulong n)
{
    if (mag_is_zero(x))
    {
        if (n == 0)
            mag_one(res);
        else
            mag_zero(res);
    }
    else if (mag_is_inf(x))
    {
        mag_inf(res);
    }
    else
    {
        mag_t t;
        mag_init(t);
        mag_one(t);
        mag_sub_lower(t, t, x);

        if (mag_is_zero(t))
        {
            mag_inf(res);
        }
        else
        {
            mag_pow_ui(res, x, n);
            mag_div(res, res, t);
        }

        mag_clear(t);
    }
}

