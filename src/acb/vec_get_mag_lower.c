/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
_acb_vec_get_mag_lower(mag_t bound, acb_srcptr vec, slong len)
{
    if (len < 1)
    {
        mag_inf(bound);
    }
    else
    {
        mag_t t;
        slong i;
        acb_get_mag_lower(bound, vec);
        mag_init(t);
        for (i = 1; i < len; i++)
        {
            acb_get_mag_lower(t, vec + i);
            mag_min(bound, bound, t);
        }
        mag_clear(t);
    }
}
