/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

void
acb_siegel_transform_cocycle_inv(acb_mat_t w, acb_mat_t c, acb_mat_t cinv,
    const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t alpha;
    fmpz_mat_t beta;
    acb_mat_t x, num;
    int r;

    fmpz_mat_window_init(alpha, mat, 0, 0, g, g);
    fmpz_mat_window_init(beta, mat, 0, g, g, 2 * g);
    acb_mat_init(x, g, g);
    acb_mat_init(num, g, g);

    acb_mat_set_fmpz_mat(x, alpha);
    acb_mat_mul(num, x, tau, prec);
    acb_mat_set_fmpz_mat(x, beta);
    acb_mat_add(num, num, x, prec);

    acb_siegel_cocycle(c, mat, tau, prec);
    r = acb_mat_inv(cinv, c, prec);
    if (!r)
    {
        acb_mat_indeterminate(cinv);
    }
    acb_mat_mul(w, num, cinv, prec);

    fmpz_mat_window_clear(alpha);
    fmpz_mat_window_clear(beta);
    acb_mat_clear(x);
    acb_mat_clear(num);
}
