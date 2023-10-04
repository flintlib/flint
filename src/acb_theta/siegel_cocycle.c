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
acb_siegel_cocycle(acb_mat_t c, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t block;
    acb_mat_t r;

    fmpz_mat_init(block, g, g);
    acb_mat_init(r, g, g);

    sp2gz_get_c(block, mat);
    acb_mat_set_fmpz_mat(c, block);
    acb_mat_mul(c, c, tau, prec);
    sp2gz_get_d(block, mat);
    acb_mat_set_fmpz_mat(r, block);
    acb_mat_add(c, c, r, prec);

    fmpz_mat_clear(block);
    acb_mat_clear(r);
}
