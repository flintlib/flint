/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_theta_g2_chi5(acb_t res, acb_srcptr th, slong prec)
{
    slong g = 2;
    slong n = 1 << (2 * g);
    ulong ab;
    acb_t t;

    acb_init(t);
    acb_one(t);

    for (ab = 0; ab < n; ab++)
    {
        if (acb_theta_char_is_even(ab, g))
        {
            acb_mul(t, t, &th[ab], prec);
        }
    }
    acb_neg(res, t);
    acb_mul_2exp_si(res, res, -6);

    acb_clear(t);
}
