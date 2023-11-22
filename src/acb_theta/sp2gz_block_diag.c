/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "acb_theta.h"

void
sp2gz_block_diag(fmpz_mat_t mat, const fmpz_mat_t U)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t D, zero;
    fmpz_t den;

    fmpz_mat_init(D, g, g);
    fmpz_mat_init(zero, g, g);
    fmpz_init(den);

    fmpz_mat_inv(D, den, U);
    fmpz_mat_transpose(D, D);

    if (!fmpz_is_one(den))
    {
        fmpz_neg(den, den);
        fmpz_mat_neg(D, D);
    }

    sp2gz_set_blocks(mat, U, zero, zero, D);

    fmpz_mat_clear(D);
    fmpz_mat_clear(zero);
    fmpz_clear(den);
}
