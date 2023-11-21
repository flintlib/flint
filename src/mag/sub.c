/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"
#include "arf.h"

void
mag_sub(mag_t z, const mag_t x, const mag_t y)
{
    if (mag_is_special(x) || mag_is_special(y))
    {
        if (mag_is_inf(x))
            mag_inf(z);
        else if (mag_is_inf(y))
            mag_zero(z);
        else
            mag_set(z, x);   /* x - 0 = x;  max(0 - y, 0) = 0 */
    }
    else
    {
        /* this should eventually have a proper implementation... */
        arf_t t, u;
        arf_init(t);
        arf_init(u);

        arf_set_mag(t, x);
        arf_set_mag(u, y);

        arf_sub(t, t, u, MAG_BITS, ARF_RND_UP);

        if (arf_sgn(t) >= 0)
            arf_get_mag(z, t);
        else
            mag_zero(z);

        arf_clear(t);
        arf_clear(u);
    }
}

