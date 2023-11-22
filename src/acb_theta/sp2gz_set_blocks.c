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
sp2gz_set_blocks(fmpz_mat_t mat, const fmpz_mat_t alpha, const fmpz_mat_t beta,
    const fmpz_mat_t gamma, const fmpz_mat_t delta)
{
    slong g = sp2gz_dim(mat);
    slong j, k;

    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            fmpz_set(fmpz_mat_entry(mat, j, k), fmpz_mat_entry(alpha, j, k));
            fmpz_set(fmpz_mat_entry(mat, j, k + g), fmpz_mat_entry(beta, j, k));
            fmpz_set(fmpz_mat_entry(mat, j + g, k), fmpz_mat_entry(gamma, j, k));
            fmpz_set(fmpz_mat_entry(mat, j + g, k + g), fmpz_mat_entry(delta, j, k));
        }
    }
}
