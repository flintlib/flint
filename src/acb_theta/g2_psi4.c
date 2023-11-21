/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_theta_g2_psi4(acb_t res, acb_srcptr th2, slong prec)
{
    slong g = 2;
    ulong ab;
    acb_t s, t;

    acb_init(s);
    acb_init(t);

    for (ab = 0; ab < (1 << (2 * g)); ab++)
    {
        if (acb_theta_char_is_even(ab, g))
        {
            acb_pow_ui(t, &th2[ab], 4, prec);
            acb_add(s, s, t, prec);
        }
    }
    acb_mul_2exp_si(res, s, -2);

    acb_clear(s);
    acb_clear(t);
}
