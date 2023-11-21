/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
sp2gz_inv(fmpz_mat_t inv, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t j;

    fmpz_mat_init(j, 2 * g, 2 * g);

    sp2gz_j(j);
    fmpz_mat_transpose(inv, mat);
    fmpz_mat_mul(inv, inv, j);
    fmpz_mat_mul(inv, j, inv);
    fmpz_mat_neg(inv, inv);

    fmpz_mat_clear(j);
}
