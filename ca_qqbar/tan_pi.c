/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_tan_pi(ca_qqbar_t res, slong p, ulong q)
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
        ca_qqbar_zero(res);
    }
    else if (q == 2)
    {
        flint_printf("ca_qqbar_tan_pi: division by zero\n");
        flint_abort();
    }
    else if (q == 4)
    {
        if (p % 4 == 1 || p % 4 == -3)
            ca_qqbar_one(res);
        else
            ca_qqbar_set_si(res, -1);
    }
    else if (q == 3)
    {
        ca_qqbar_set_ui(res, 3);
        ca_qqbar_sqrt(res, res);
        if (p % 3 == -1 || p % 3 == 2)
            ca_qqbar_neg(res, res);
    }
    else if (q == 6)
    {
        ca_qqbar_set_ui(res, 3);
        ca_qqbar_sqrt(res, res);
        ca_qqbar_inv(res, res);
        if (p % 6 == -1 || p % 6 == 5)
            ca_qqbar_neg(res, res);
    }
    else
    {
        ca_qqbar_t t;
        ca_qqbar_init(t);

        ca_qqbar_exp_pi_i(res, 2 * p, q);
        ca_qqbar_add_ui(res, res, 1);
        ca_qqbar_inv(res, res);
        ca_qqbar_mul_2exp_si(res, res, 1);
        ca_qqbar_sub_ui(res, res, 1);
        ca_qqbar_i(t);
        ca_qqbar_mul(res, res, t);

        arb_zero(acb_imagref(CA_QQBAR_ENCLOSURE(res)));

        ca_qqbar_clear(t);
    }
}

