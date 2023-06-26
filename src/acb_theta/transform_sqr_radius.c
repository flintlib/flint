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
acb_theta_transform_sqr_radius(arf_t rho, const arf_t r, acb_srcptr th2,
                               const fmpz_mat_t mat, slong prec)
{
    ulong ab_0, ab;
    fmpz_t eps;
    arb_t abs_0, abs, t;
    arf_t bound, max, res;
    slong g = fmpz_mat_nrows(mat) / 2;
    slong n = 1 << g;
    slong k;

    fmpz_init(eps);
    arb_init(abs_0);
    arb_init(abs);
    arb_init(t);
    arf_init(bound);
    arf_init(max);
    arf_init(res);

    ab_0 = acb_theta_transform_image_char(eps, 0, mat);

    /* Compute suitable radius for duplicated values */
    acb_abs(abs_0, &th2[ab_0], prec);
    arf_pos_inf(res);

    for (k = 1; k < n; k++)
    {
        ab = acb_theta_transform_image_char(eps, k, mat);
        acb_abs(abs, &th2[ab], prec);

        arb_one(t);
        arb_add_arf(t, t, r, prec);
        arb_mul(t, t, abs_0, prec);
        arb_add(t, t, abs, prec);
        arb_div(t, abs_0, t, prec);
        arb_mul(t, t, abs_0, prec);
        arb_mul_arf(t, t, r, prec);

        arb_get_lbound_arf(bound, t, prec);
        arf_min(res, res, bound);
    }

    arf_set(rho, res);

    fmpz_clear(eps);
    arb_clear(abs_0);
    arb_clear(abs);
    arb_clear(t);
    arf_clear(bound);
    arf_clear(max);
    arf_clear(res);
}
