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
fmpz_mat_set_abcd(fmpz_mat_t mat, const fmpz_mat_t a, const fmpz_mat_t b,
                  const fmpz_mat_t c, const fmpz_mat_t d)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    slong j, k;

    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            fmpz_set(fmpz_mat_entry(mat, j, k), fmpz_mat_entry(a, j, k));
            fmpz_set(fmpz_mat_entry(mat, j, k + g), fmpz_mat_entry(b, j, k));
            fmpz_set(fmpz_mat_entry(mat, j + g, k), fmpz_mat_entry(c, j, k));
            fmpz_set(fmpz_mat_entry(mat, j + g, k + g),
                     fmpz_mat_entry(d, j, k));
        }
    }
}
