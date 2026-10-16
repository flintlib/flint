/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_mat.h"

slong
nmod_mat_left_nullspace(nmod_mat_t X, const nmod_mat_t A)
{
    slong m = A->r;
    slong n = A->c;
    slong nullity, i, j;
    nmod_mat_t At, Y;

    /*
        Left nullspace of A equals right nullspace of A^T, expressed
        as row vectors.  We defer to nmod_mat_nullspace on A^T (which
        returns column vectors) and transpose the result.
    */
    nmod_mat_init(At, n, m, A->mod.n);
    nmod_mat_transpose(At, A);

    nmod_mat_init(Y, m, m, A->mod.n);
    nullity = nmod_mat_nullspace(Y, At);

    nmod_mat_zero(X);
    for (i = 0; i < nullity; i++)
        for (j = 0; j < m; j++)
            nmod_mat_entry(X, i, j) = nmod_mat_entry(Y, j, i);

    nmod_mat_clear(Y);
    nmod_mat_clear(At);

    return nullity;
}
