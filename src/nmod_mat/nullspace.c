/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

slong
nmod_mat_nullspace(nmod_mat_t X, const nmod_mat_t A)
{
    slong m = A->r, n = A->c, rank, nullity;
    slong * P, * Q;
    nmod_mat_t LU, Y;

    P = flint_malloc((m + 1) * sizeof(slong));
    Q = flint_malloc((n + 1) * sizeof(slong));
    nmod_mat_init_set(LU, A);

    rank = nmod_mat_ple(P, Q, LU, 0);
    nullity = n - rank;

    /* X may be aliased with A, which is no longer read; its first nullity
       columns are written by nmod_mat_ple_nullspace, the others zeroed */
    if (nullity != 0)
    {
        nmod_mat_window_init(Y, X, 0, 0, n, nullity);
        nmod_mat_ple_nullspace(Y, LU, Q, rank);
        nmod_mat_window_clear(Y);
    }
    if (X->c > nullity)
    {
        nmod_mat_window_init(Y, X, 0, nullity, X->r, X->c);
        nmod_mat_zero(Y);
        nmod_mat_window_clear(Y);
    }

    nmod_mat_clear(LU);
    flint_free(P);
    flint_free(Q);

    return nullity;
}
