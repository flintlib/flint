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
sp2gz_get_a(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    slong j, k;

    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            fmpz_set(fmpz_mat_entry(res, j, k), fmpz_mat_entry(mat, j, k));
        }
    }
}
