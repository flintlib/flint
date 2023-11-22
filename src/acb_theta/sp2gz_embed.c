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
sp2gz_embed(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong j, k, u, v;
    slong g = sp2gz_dim(res);
    slong g1 = sp2gz_dim(mat);

    fmpz_mat_one(res);
    for (j = 0; j < 2 * g1; j++)
    {
        for (k = 0; k < 2 * g1; k++)
        {
            u = j + (j >= g1 ? g - g1 : 0);
            v = k + (k >= g1 ? g - g1 : 0);
            fmpz_set(fmpz_mat_entry(res, u, v), fmpz_mat_entry(mat, j, k));
        }
    }
}
