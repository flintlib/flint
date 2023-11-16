/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_bound_factor(mag_t C,
    acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, ulong n)
{
    mag_t t, u;
    acb_t w;
    slong i;

    if (p > q)
    {
        mag_inf(C);
        return;
    }

    mag_init(t);
    mag_init(u);
    acb_init(w);

    acb_get_mag(C, z);

    for (i = 0; i < q; i++)
    {
        acb_add_ui(w, b + i, n, MAG_BITS);

        if (arb_is_positive(acb_realref(w)))
        {
            acb_get_mag_lower(u, w);

            if (i < p)
            {
                acb_sub(w, a + i, b + i, MAG_BITS);
                acb_get_mag(t, w);
                mag_div(t, t, u);
                mag_one(u);
                mag_add(t, t, u);
                mag_mul(C, C, t);
            }
            else
            {
                mag_div(C, C, u);
            }
        }
        else
        {
            mag_inf(C);
        }
    }

    mag_one(t);
    mag_sub_lower(u, t, C);

    if (mag_is_zero(u))
        mag_inf(C);
    else
        mag_div(C, t, u);

    mag_clear(t);
    mag_clear(u);
    acb_clear(w);
}

