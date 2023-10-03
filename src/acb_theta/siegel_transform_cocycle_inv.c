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
acb_siegel_transform_cocycle_inv(acb_mat_t w, acb_mat_t c, acb_mat_t cinv,
    const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t a;
    acb_mat_t x, num;
    int r;

    fmpz_mat_init(a, g, g);
    acb_mat_init(x, g, g);
    acb_mat_init(num, g, g);

    sp2gz_get_a(a, mat);
    acb_mat_set_fmpz_mat(x, a);
    acb_mat_mul(num, x, tau, prec);
    sp2gz_get_b(a, mat);
    acb_mat_set_fmpz_mat(x, a);
    acb_mat_add(num, num, x, prec);

    acb_siegel_cocycle(c, mat, tau, prec);
    r = acb_mat_inv(cinv, c, prec);
    if (!r)
    {
        acb_mat_indeterminate(cinv);
    }
    acb_mat_mul(w, num, cinv, prec);

    fmpz_mat_clear(a);
    acb_mat_clear(x);
    acb_mat_clear(num);
}
