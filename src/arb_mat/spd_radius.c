/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

/* This uses the formula dA = A Phi(A^{-1} dP A^{-T}) where A is
   lower-triangular, P = A A^T, d denotes some scalar derivative, and
   Phi_{i,j}(M) is M_{i,j} for i>j, 1/2 M_{i,i} for i=j, and 0 for i<j */

void arb_mat_spd_radius(arf_t r, const arb_mat_t A, slong prec)
{
    slong g = arb_mat_nrows(A);
    arb_mat_t cho;
    mag_t b, binv;
    arf_t t;
    int res;
    slong k, j;

    arb_mat_init(cho, g, g);
    mag_init(b);
    mag_init(binv);
    arf_init(t);

    res = arb_mat_cho(cho, A, prec);
    if (!res)
    {
        arf_nan(r);
    }
    else
    {
        /* Get accepted deformation on cho for infinity norm */
        arf_pos_inf(r);
        for (k = 0; k < g; k++)
        {
            arb_get_lbound_arf(t, arb_mat_entry(cho, k, k), prec);
            arf_min(r, r, t);
        }
        arf_mul_2exp_si(r, r, -1);

        /* Add error to cho */
        for (k = 0; k < g; k++)
        {
            for (j = 0; j <= k; j++)
            {
                arb_add_error_arf(arb_mat_entry(cho, k, j), r);
            }
        }

        /* Get induced norms */
        arb_mat_bound_inf_norm(b, cho);
        arb_mat_inv(cho, cho, prec);
        arb_mat_bound_inf_norm(binv, cho);

        /* Propagate bound */
        arf_set_mag(t, b);
        arf_div(r, r, t, prec, ARF_RND_FLOOR);
        arf_set_mag(t, binv);
        arf_div(r, r, t, prec, ARF_RND_FLOOR);
        arf_div(r, r, t, prec, ARF_RND_FLOOR);
    }

    arb_mat_clear(cho);
    mag_clear(b);
    mag_clear(binv);
    arf_clear(t);
}
