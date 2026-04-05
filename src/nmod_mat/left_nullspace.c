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
    slong m = nmod_mat_nrows(A);
    slong n = nmod_mat_ncols(A);
    slong rank_A, nullity, i, j;
    nmod_mat_t aug;

    /*
        Augment [A | I_m] and row-reduce.  In RREF, pivots in
        columns 0..n-1 correspond to rank(A); the remaining rows
        have their A part entirely zero and the I part gives the
        left nullspace vectors.
    */
    nmod_mat_init(aug, m, n + m, A->mod.n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            nmod_mat_entry(aug, i, j) = nmod_mat_entry(A, i, j);
        nmod_mat_entry(aug, i, n + i) = UWORD(1);
    }

    nmod_mat_rref(aug);

    /*
        RREF pivot columns are in increasing order, so rows with
        pivots in A columns (0..n-1) come first.  Find rank(A) by
        locating the first row whose A part is entirely zero.
    */
    for (rank_A = 0; rank_A < m; rank_A++)
    {
        int found = 0;
        for (j = 0; j < n; j++)
        {
            if (nmod_mat_entry(aug, rank_A, j) != UWORD(0))
            {
                found = 1;
                break;
            }
        }
        if (!found)
            break;
    }

    nullity = m - rank_A;

    nmod_mat_zero(X);

    for (j = 0; j < nullity; j++)
    {
        for (i = 0; i < m; i++)
            nmod_mat_entry(X, i, j) =
                nmod_mat_entry(aug, rank_A + j, n + i);
    }

    nmod_mat_clear(aug);

    return nullity;
}
