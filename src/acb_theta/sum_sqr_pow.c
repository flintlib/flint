/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_sum_sqr_pow(acb_ptr * sqr_pow, const acb_mat_t exp_tau, const acb_theta_eld_t E, slong prec)
{
    slong g = acb_mat_nrows(exp_tau);
    acb_t c, dc, ddc;
    slong k, j;

    acb_init(c);
    acb_init(dc);
    acb_init(ddc);

    /* Addition chains do not make a huge difference here. */
    for (k = 0; k < g; k++)
    {
        acb_one(c);
        acb_set(dc, acb_mat_entry(exp_tau, k, k));
        acb_sqr(ddc, dc, prec);
        for (j = 0; j <= acb_theta_eld_box(E, k); j++)
        {
            acb_set(&sqr_pow[k][j], c);
            acb_mul(c, c, dc, prec);
            acb_mul(dc, dc, ddc, prec);
        }
    }

    acb_clear(c);
    acb_clear(dc);
    acb_clear(ddc);
}
