/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int acb_siegel_is_reduced(const acb_mat_t tau, slong tol_exp, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t im;
    arb_t abs, t, u;
    slong j, k;
    int res = 1;

    arb_mat_init(im, g, g);
    arb_init(abs);
    arb_init(t);
    arb_init(u);

    arb_one(u);
    arb_mul_2exp_si(u, u, tol_exp);

    arb_set_si(t, 3);
    arb_sqrt(t, t, prec);
    arb_mul_2exp_si(t, t, -1);
    arb_sub(t, t, u, prec);
    if (!arb_gt(acb_imagref(acb_mat_entry(tau, 0, 0)), t))
    {
        res = 0;
    }

    arb_one(t);
    arb_mul_2exp_si(t, t, -1);
    arb_add(t, t, u, prec);
    for (j = 0; (j < g) && res; j++)
    {
        for (k = 0; (k < g) && res; k++)
        {
            arb_abs(abs, acb_realref(acb_mat_entry(tau, j, k)));
            if (!arb_lt(abs, t))
            {
                res = 0;
            }
        }
    }

    if (res)
    {
        acb_mat_get_imag(im, tau);
        res = arb_mat_spd_is_lll_reduced(im, tol_exp, prec);
    }

    arb_mat_clear(im);
    arb_clear(abs);
    arb_clear(t);
    arb_clear(u);
    return res;
}
