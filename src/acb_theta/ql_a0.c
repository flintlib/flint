/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static slong
acb_theta_ql_split(const arb_mat_t cho)
{
    slong g = arb_mat_nrows(cho);
    arb_t cmp;
    slong k;

    arb_init(cmp);

    for (k = g - 1; k >= 1; k--)
    {
        arb_mul_2exp_si(cmp, arb_mat_entry(cho, k - 1, k - 1),
            ACB_THETA_QL_SPLIT);
        if (arb_lt(cmp, arb_mat_entry(cho, k, k)))
        {
            break;
        }
    }

    arb_clear(cmp);
    return k;
}

int
acb_theta_ql_a0(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t cho;
    slong split, nb_steps;
    int res;

    arb_mat_init(cho, g, g);

    acb_theta_eld_cho(cho, tau, ACB_THETA_LOW_PREC);
    split = acb_theta_ql_split(cho);
    nb_steps = acb_theta_ql_nb_steps(cho, split, prec);

    res = acb_theta_ql_a0_steps(r, t, z, dist0, dist, tau, nb_steps, split,
        guard, prec, &acb_theta_ql_a0);

    arb_mat_clear(cho);
    return res;
}
