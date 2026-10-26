/*
    Copyright (C) 2026 Edgar Costa
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

slong
nmod_mat_left_nullspace(nmod_mat_t X, const nmod_mat_t A)
{
    slong m = A->r, rank, nullity;
    slong * P;
    nmod_mat_t LU, Y;

    P = flint_malloc((m + 1) * sizeof(slong));
    nmod_mat_init_set(LU, A);

    /* the row order of nmod_mat_ple reveals the row rank profile, so that
       no transposition is needed */
    rank = nmod_mat_ple(P, NULL, LU, 0);
    nullity = m - rank;

    /* the first nullity rows and m columns of X are written by
       nmod_mat_ple_left_nullspace, the others zeroed */
    if (nullity != 0)
    {
        nmod_mat_window_init(Y, X, 0, 0, nullity, m);
        nmod_mat_ple_left_nullspace(Y, LU, P, rank);
        nmod_mat_window_clear(Y);
    }
    if (X->c > m && nullity != 0)
    {
        nmod_mat_window_init(Y, X, 0, m, nullity, X->c);
        nmod_mat_zero(Y);
        nmod_mat_window_clear(Y);
    }
    if (X->r > nullity)
    {
        nmod_mat_window_init(Y, X, nullity, 0, X->r, X->c);
        nmod_mat_zero(Y);
        nmod_mat_window_clear(Y);
    }

    nmod_mat_clear(LU);
    flint_free(P);

    return nullity;
}
