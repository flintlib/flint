/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_get_mag_lower(mag_t z, const acb_t x)
{
    if (arb_is_zero(acb_imagref(x)))
    {
        arb_get_mag_lower(z, acb_realref(x));
    }
    else if (arb_is_zero(acb_realref(x)))
    {
        arb_get_mag_lower(z, acb_imagref(x));
    }
    else
    {
        mag_t t;
        mag_init(t);

        arb_get_mag_lower(t, acb_realref(x));
        arb_get_mag_lower(z, acb_imagref(x));

        mag_mul_lower(t, t, t);
        mag_mul_lower(z, z, z);
        mag_add_lower(z, z, t);
        mag_sqrt_lower(z, z);

        mag_clear(t);
    }
}

