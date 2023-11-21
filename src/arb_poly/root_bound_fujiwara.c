/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_root_bound_fujiwara(mag_t bound, arb_srcptr poly, slong len)
{
    mag_t t, u, v;
    slong i;

    if (len <= 1)
    {
        mag_inf(bound);
        return;
    }

    mag_init(t);
    mag_init(u);
    mag_init(v);

    /* u = 1/leading */
    arb_get_mag_lower(t, poly + len - 1);
    mag_one(u);
    mag_div(u, u, t);

    mag_zero(v);

    for (i = 0; i < len - 1; i++)
    {
        arb_get_mag(t, poly + len - 2 - i);
        mag_mul(t, t, u);

        if (i == len - 2)
            mag_mul_2exp_si(t, t, -1);

        mag_root(t, t, i + 1);
        mag_max(v, v, t);
    }

    mag_mul_2exp_si(bound, v, 1);

    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
}

void
arb_poly_root_bound_fujiwara(mag_t bound, arb_poly_t poly)
{
    _arb_poly_root_bound_fujiwara(bound, poly->coeffs, poly->length);
}

