/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


static __inline__ int
nmod_mat_pivot(nmod_mat_t A, slong * P, slong start_row, slong col)
{
    slong j, t;
    mp_ptr u;

    if (nmod_mat_entry(A, start_row, col) != 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (nmod_mat_entry(A, j, col) != 0)
        {
            u = A->rows[j];
            A->rows[j] = A->rows[start_row];
            A->rows[start_row] = u;

            t = P[j];
            P[j] = P[start_row];
            P[start_row] = t;

            return -1;
        }
    }
    return 0;
}


slong 
nmod_mat_lu_classical(slong * P, nmod_mat_t A, int rank_check)
{
    mp_limb_t d, e, **a;
    nmod_t mod;
    slong i, m, n, rank, length, row, col;

    m = A->r;
    n = A->c;
    a = A->rows;
    mod = A->mod;

    rank = row = col = 0;

    for (i = 0; i < m; i++)
        P[i] = i;

    while (row < m && col < n)
    {
        if (nmod_mat_pivot(A, P, row, col) == 0)
        {
            if (rank_check)
                return 0;
            col++;
            continue;
        }

        rank++;

        d = a[row][col];
        d = n_invmod(d, mod.n);
        length = n - col - 1;

        for (i = row + 1; i < m; i++)
        {
            e = n_mulmod2_preinv(a[i][col], d, mod.n, mod.ninv);
            if (length != 0)
                _nmod_vec_scalar_addmul_nmod(a[i] + col + 1,
                    a[row] + col + 1, length, nmod_neg(e, mod), mod);

            a[i][col] = 0;
            a[i][rank - 1] = e;
        }
        row++;
        col++;
    }

    return rank;
}
