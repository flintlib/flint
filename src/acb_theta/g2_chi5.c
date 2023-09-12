/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_g2_chi5(acb_t r, acb_srcptr th, slong prec)
{
    slong g = 2;
    slong n = 1 << (2 * g);
    ulong ab;
    acb_t res;

    acb_init(res);
    acb_one(res);

    for (ab = 0; ab < n; ab++)
    {
        if (theta_char_is_even(ab, g))
        {
            acb_mul(res, res, &th[ab], prec);
        }
    }
    acb_neg(r, res);
    acb_scalar_mul_2exp_si(r, r, -6);

    acb_clear(res);
}
