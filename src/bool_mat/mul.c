/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

void
bool_mat_mul(bool_mat_t C, const bool_mat_t A, const bool_mat_t B)
{
    slong ar, ac, br, bc, i, j, k;

    ar = bool_mat_nrows(A);
    ac = bool_mat_ncols(A);
    br = bool_mat_nrows(B);
    bc = bool_mat_ncols(B);

    if (ac != br || ar != bool_mat_nrows(C) || bc != bool_mat_ncols(C))
    {
        flint_throw(FLINT_ERROR, "bool_mat_mul: incompatible dimensions\n");
    }

    if (br == 0)
    {
        bool_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        bool_mat_t T;
        bool_mat_init(T, ar, bc);
        bool_mat_mul(T, A, B);
        bool_mat_swap(T, C);
        bool_mat_clear(T);
        return;
    }

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            int any = 0;
            for (k = 0; k < br && !any; k++)
                any |= (bool_mat_get_entry(A, i, k) &
                        bool_mat_get_entry(B, k, j));
            bool_mat_set_entry(C, i, j, any);
        }
    }
}
