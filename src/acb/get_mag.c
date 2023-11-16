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
acb_get_mag(mag_t u, const acb_t z)
{
    if (arb_is_zero(acb_imagref(z)))
    {
        arb_get_mag(u, acb_realref(z));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        arb_get_mag(u, acb_imagref(z));
    }
    else
    {
        mag_t v;
        mag_init(v);

        arb_get_mag(u, acb_realref(z));
        arb_get_mag(v, acb_imagref(z));

        mag_mul(u, u, u);
        mag_addmul(u, v, v);
        mag_sqrt(u, u);

        mag_clear(v);
    }
}

