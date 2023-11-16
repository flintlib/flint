/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_hurwitz_zeta_uiui(mag_t res, ulong s, ulong a)
{
    if (s <= 1 || a == 0)
    {
        mag_inf(res);
    }
    else
    {
        mag_t t, u;
        mag_init(t);
        mag_init(u);

        mag_one(t);
        mag_set_ui_lower(u, a);
        mag_pow_ui_lower(u, u, s - 1);
        mag_mul_ui_lower(res, u, a);
        mag_div(res, t, res);
        mag_mul_ui_lower(u, u, s - 1);
        mag_div(u, t, u);
        mag_add(res, res, u);

        mag_clear(t);
        mag_clear(u);
    }
}

