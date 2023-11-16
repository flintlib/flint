/*
    Copyright 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_eig_global_enclosure(mag_t eps, const acb_mat_t A, acb_srcptr E, const acb_mat_t R, slong prec)
{
    acb_mat_t Y, R1, R2;
    slong i, j, n;
    mag_t r1, r2;

    n = acb_mat_nrows(A);

    acb_mat_init(Y, n, n);
    acb_mat_init(R1, n, n);
    acb_mat_init(R2, n, n);
    mag_init(r1);
    mag_init(r2);

    /* Y ~= inv(R) */
    acb_mat_one(R1);
    acb_mat_approx_solve(Y, R, R1, prec);

    /* R2 = Y*R - I */
    acb_mat_mul(R2, Y, R, prec);
    for (i = 0; i < n; i++)
        acb_sub_ui(acb_mat_entry(R2, i, i), acb_mat_entry(R2, i, i), 1, prec);

    acb_mat_bound_inf_norm(r2, R2);

    if (mag_cmp_2exp_si(r2, 0) < 0)
    {
        /* R1 = Y*(AR - RD) */
        acb_mat_mul(R2, A, R, prec);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                acb_submul(acb_mat_entry(R2, i, j), acb_mat_entry(R, i, j), E + j, prec);
        acb_mat_mul(R1, Y, R2, prec);

        acb_mat_bound_inf_norm(r1, R1);
        mag_geom_series(r2, r2, 0);
        mag_mul(eps, r1, r2);
    }
    else
    {
        mag_inf(eps);
    }

    acb_mat_clear(R1);
    acb_mat_clear(R2);
    acb_mat_clear(Y);
    mag_clear(r1);
    mag_clear(r2);
}
