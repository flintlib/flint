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
acb_theta_g2_psi4(acb_t r, acb_srcptr th2, slong prec)
{
    slong g = 2;
    ulong ab;
    acb_t res, aux;

    acb_init(res);
    acb_init(aux);

    for (ab = 0; ab < (1 << (2 * g)); ab++)
    {
        if (acb_theta_char_is_even(ab, g))
        {
            acb_pow_ui(aux, &th2[ab], 4, prec);
            acb_add(res, res, aux, prec);
        }
    }
    acb_mul_2exp_si(r, res, -2);

    acb_clear(res);
    acb_clear(aux);
}
