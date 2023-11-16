/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "arb_mat.h"

int
arb_mat_spd_get_fmpz_mat(fmpz_mat_t B, const arb_mat_t A, slong prec)
{
    slong j, k;
    slong g = arb_mat_nrows(A);
    arb_t z;
    int res = 1;

    arb_init(z);

    for (j = 0; (j < g) && res; j++)
    {
        for (k = j; (k < g) && res; k++)
        {
            res = arb_intersection(z, arb_mat_entry(A, j, k),
                arb_mat_entry(A, k, j), prec);
            arf_get_fmpz_fixed_si(fmpz_mat_entry(B, j, k), arb_midref(z), -prec);
            fmpz_set(fmpz_mat_entry(B, k, j), fmpz_mat_entry(B, j, k));
        }
    }
    res = res && fmpz_mat_is_spd(B);

    arb_clear(z);
    return res;
}

