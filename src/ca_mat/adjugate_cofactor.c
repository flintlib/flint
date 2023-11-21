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
ca_mat_adjugate_cofactor(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_mat_t T;
    slong i, j, n, a, b;
    ca_t t, zero;

    n = ca_mat_nrows(A);

    if (n == 0)
    {
        ca_one(det, ctx);
        return;
    }

    if (n == 1)
    {
        ca_set(det, ca_mat_entry(A, 0, 0), ctx);
        ca_one(ca_mat_entry(adj, 0, 0), ctx);
        return;
    }

    if (n == 2)
    {
        ca_t t, u;
        ca_init(t, ctx);
        ca_init(u, ctx);
        ca_mul(t, ca_mat_entry(A, 0, 0), ca_mat_entry(A, 1, 1), ctx);
        ca_mul(u, ca_mat_entry(A, 0, 1), ca_mat_entry(A, 1, 0), ctx);
        ca_set(ca_mat_entry(adj, 0, 0), ca_mat_entry(A, 1, 1), ctx);
        ca_neg(ca_mat_entry(adj, 0, 1), ca_mat_entry(A, 0, 1), ctx);
        ca_neg(ca_mat_entry(adj, 1, 0), ca_mat_entry(A, 1, 0), ctx);
        ca_set(ca_mat_entry(adj, 1, 1), ca_mat_entry(A, 0, 0), ctx);
        ca_sub(det, t, u, ctx);
        ca_clear(t, ctx);
        ca_clear(u, ctx);
        return;
    }

    if (adj == A)
    {
        ca_mat_init(T, n, n, ctx);
        ca_mat_adjugate_cofactor(T, det, A, ctx);
        ca_mat_swap(adj, T, ctx);
        ca_mat_clear(T, ctx);
        return;
    }

    ca_mat_init(T, n - 1, n - 1, ctx);
    ca_init(zero, ctx);
    ca_init(t, ctx);

    ca_zero(det, ctx);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (a = 0; a < n; a++)
            {
                for (b = 0; b < n; b++)
                {
                    if (a != i && b != j)
                    {
                        *ca_mat_entry(T, a - (a > i), b - (b > j)) = *ca_mat_entry(A, a, b);
                    }
                }
            }

            ca_mat_det(ca_mat_entry(adj, i, j), T, ctx);
            if ((i + j) & 1)
                ca_neg(ca_mat_entry(adj, i, j), ca_mat_entry(adj, i, j), ctx);

            if (i == 0)
            {
                ca_mul(t, ca_mat_entry(adj, i, j), ca_mat_entry(A, i, j), ctx);
                ca_add(det, det, t, ctx);
            }
        }
    }

    ca_mat_transpose(adj, adj, ctx);

    for (i = 0; i < n - 1; i++)
        for (j = 0; j < n - 1; j++)
            *ca_mat_entry(T, i, j) = *zero;
    ca_mat_clear(T, ctx);
    ca_clear(t, ctx);
}
