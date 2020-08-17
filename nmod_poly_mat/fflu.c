/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

#define E(j,k) nmod_poly_mat_entry(B,j,k)

static __inline__ void
nmod_poly_mat_swap_rows(nmod_poly_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        nmod_poly_struct * u;
        slong t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u; 
    }
}

slong
nmod_poly_mat_fflu(nmod_poly_mat_t B, nmod_poly_t den, slong * perm,
    const nmod_poly_mat_t A, int rank_check)
{
    nmod_poly_t t;
    slong m, n, j, k, rank, r, pivot_row, pivot_col;

    if (nmod_poly_mat_is_empty(A))
    {
        nmod_poly_one(den);
        return 0;
    }

    nmod_poly_mat_set(B, A);
    m = B->r;
    n = B->c;
    rank = pivot_row = pivot_col = 0;

    nmod_poly_init(t, nmod_poly_mat_modulus(A));

    while (pivot_row < m && pivot_col < n)
    {
        r = nmod_poly_mat_find_pivot_partial(B, pivot_row, m, pivot_col);

        if (r == -1)
        {
            if (rank_check)
            {
                nmod_poly_zero(den);
                rank = 0;
                break;
            }
            pivot_col++;
            continue;
        }
        else if (r != pivot_row)
            nmod_poly_mat_swap_rows(B, perm, pivot_row, r);

        rank++;

        for (j = pivot_row + 1; j < m; j++)
        {
            for (k = pivot_col + 1; k < n; k++)
            {
                nmod_poly_mul(E(j, k), E(j, k), E(pivot_row, pivot_col));
                nmod_poly_mul(t, E(j, pivot_col), E(pivot_row, k));
                nmod_poly_sub(E(j, k), E(j, k), t);
                if (pivot_row > 0)
                    nmod_poly_div(E(j, k), E(j, k), den);
            }
        }

        nmod_poly_set(den, E(pivot_row, pivot_col));
        pivot_row++;
        pivot_col++;
    }

    nmod_poly_clear(t);
    return rank;
}
