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
acb_siegel_cocycle(acb_mat_t c, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t gamma, delta;
    acb_mat_t r;

    fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
    fmpz_mat_window_init(delta, mat, g, g, 2 * g, 2 * g);
    acb_mat_init(r, g, g);

    acb_mat_set_fmpz_mat(c, gamma);
    acb_mat_set_fmpz_mat(r, delta);
    acb_mat_mul(c, c, tau, prec);
    acb_mat_add(c, c, r, prec);

    fmpz_mat_window_clear(gamma);
    fmpz_mat_window_clear(delta);
    acb_mat_clear(r);
}
