/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
sp2gz_is_gsp(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t r, J;
    int res;

    fmpz_mat_init(r, 2 * g, 2 * g);
    fmpz_mat_init(J, 2 * g, 2 * g);
    sp2gz_j(J);

    fmpz_mat_transpose(r, mat);
    fmpz_mat_mul(r, r, J);
    fmpz_mat_mul(r, r, mat);
    fmpz_mat_mul(r, r, J);

    res = sp2gz_is_scalar(r) && !fmpz_is_zero(fmpz_mat_entry(r, 0, 0));

    fmpz_mat_clear(r);
    fmpz_mat_clear(J);
    return res;
}
