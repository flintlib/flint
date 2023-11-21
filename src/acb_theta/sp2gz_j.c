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
sp2gz_j(fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t zero, one, minus_one;

    fmpz_mat_init(zero, g, g);
    fmpz_mat_init(one, g, g);
    fmpz_mat_init(minus_one, g, g);

    fmpz_mat_one(one);
    fmpz_mat_neg(minus_one, one);
    sp2gz_set_blocks(mat, zero, one, minus_one, zero);

    fmpz_mat_clear(zero);
    fmpz_mat_clear(one);
    fmpz_mat_clear(minus_one);
}
