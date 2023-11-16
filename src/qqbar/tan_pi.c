/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "qqbar.h"

int
qqbar_tan_pi(qqbar_t res, slong p, ulong q)
{
    slong g;

    g = n_gcd(FLINT_ABS(p), q);

    if (g != 1)
    {
        p /= g;
        q /= g;
    }

    if (q == 1)
    {
        qqbar_zero(res);
    }
    else if (q == 2)
    {
        return 0;
    }
    else if (q == 4)
    {
        if (p % 4 == 1 || p % 4 == -3)
            qqbar_one(res);
        else
            qqbar_set_si(res, -1);
    }
    else if (q == 3)
    {
        qqbar_set_ui(res, 3);
        qqbar_sqrt(res, res);
        if (p % 3 == -1 || p % 3 == 2)
            qqbar_neg(res, res);
    }
    else if (q == 6)
    {
        qqbar_set_ui(res, 3);
        qqbar_sqrt(res, res);
        qqbar_inv(res, res);
        if (p % 6 == -1 || p % 6 == 5)
            qqbar_neg(res, res);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);

        qqbar_exp_pi_i(res, 2 * p, q);
        qqbar_add_ui(res, res, 1);
        qqbar_inv(res, res);
        qqbar_mul_2exp_si(res, res, 1);
        qqbar_sub_ui(res, res, 1);
        qqbar_i(t);
        qqbar_mul(res, res, t);

        arb_zero(acb_imagref(QQBAR_ENCLOSURE(res)));

        qqbar_clear(t);
    }

    return 1;
}
