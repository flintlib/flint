/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_mul_classical(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    slong ar, ac, br, bc, i, j, k;
    ca_t t;

    ar = ca_mat_nrows(A);
    ac = ca_mat_ncols(A);
    br = ca_mat_nrows(B);
    bc = ca_mat_ncols(B);

    if (ac != br || ar != ca_mat_nrows(C) || bc != ca_mat_ncols(C))
    {
        flint_throw(FLINT_ERROR, "ca_mat_mul_classical: incompatible dimensions\n");
    }

    if (br == 0)
    {
        ca_mat_zero(C, ctx);
        return;
    }

    if (A == C || B == C)
    {
        ca_mat_t T;
        ca_mat_init(T, ar, bc, ctx);
        ca_mat_mul(T, A, B, ctx);
        ca_mat_swap(T, C, ctx);
        ca_mat_clear(T, ctx);
        return;
    }

    ca_init(t, ctx);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            ca_mul(ca_mat_entry(C, i, j),
                      ca_mat_entry(A, i, 0),
                      ca_mat_entry(B, 0, j), ctx);

            for (k = 1; k < br; k++)
            {
                ca_mul(t, ca_mat_entry(A, i, k), ca_mat_entry(B, k, j), ctx);
                ca_add(ca_mat_entry(C, i, j), ca_mat_entry(C, i, j), t, ctx);
            }
        }
    }

    ca_clear(t, ctx);
}
