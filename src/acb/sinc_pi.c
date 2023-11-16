/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_sinc_pi(acb_t res, const acb_t x, slong prec)
{
    mag_t m;
    acb_t t;

    if (acb_is_zero(x))
    {
        acb_one(res);
        return;
    }

    mag_init(m);
    acb_init(t);

    acb_get_mag_lower(m, x);

    if (mag_cmp_2exp_si(m, -1) > 0)
    {
        acb_const_pi(t, prec + 4);
        acb_mul(t, t, x, prec + 4);
        acb_sin_pi(res, x, prec + 4);
        acb_div(res, res, t, prec);
    }
    else
    {
        acb_const_pi(t, prec + 4);
        acb_mul(t, t, x, prec + 4);
        acb_sinc(res, t, prec);
    }

    mag_clear(m);
    acb_clear(t);
}

