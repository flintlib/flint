/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

slong fmpz_sparse_mat_hnf_xgcd(fmpz_sparse_mat_t M)
{
    slong i, r, pr, pc, rank, remr, nnz;
    slong *P; /* Row permutation */
    slong *irows; /* Set of rows incident on a given column */
    slong nnp; /* Number of rows which are not pivots */
    slong npp; /* Number of rows which are previous pivots */
    fmpz_sparse_mat_with_transpose_t MT;
    hashmap_struct *hcol; /* Virtual L^t and one of its rows */
    fmpz_t g, a, b, one;

    if (M->r == 0 || M->c == 0) return 0;
    fmpz_init(g);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init_set_ui(one, UWORD(1));

    /* Construct virtual transpose */
    _fmpz_sparse_mat_with_transpose_init(MT, M);
    
    /* Set up permutation */
    P = flint_malloc(M->r*sizeof(*P));
    remr = M->r;
    for (r = 0; r < M->r; ++r) 
    {
        if (!M->rows[r].nnz) P[r] = --remr; 
        else P[r] = -1;
    }
    irows = flint_malloc(M->r * sizeof(*irows));
    
    for (rank = pc = 0; pc < M->c; ++pc)
    {
        hcol = &MT->cols[pc]; nnz = hcol->num;
        if (!nnz) continue;
        pr = -1, nnp = 0, npp = 0;

        for (i = 0; i < nnz; ++i)
        {
            r = hcol->keys[i];
            if (P[r] >= 0) /* Put previous pivot rows at end of irows */
                irows[hcol->num - (++npp)] = r;
            else if (pr >= 0 && fmpz_cmpabs(LT(M, r).val, LT(M, pr).val) >= 0) 
                irows[nnp++] = r; /* Put non-pivot rows with larger entries at beginning of irows */
            else /* Replace the pivot row with the current row */
            {
                if (pr >= 0) irows[nnp++] = pr;
                pr = r;
            }
        }
        if (pr == -1) continue;

        /* Eliminate larger, non-pivot rows */
        for (i = 0; i < nnp; ++i)
        {
            r = irows[i];
            _fmpz_sparse_mat_with_transpose_gauss_elim_ext(MT, pr, r);
            if (M->rows[r].nnz == 0) P[r] = --remr;
        }

        if (fmpz_sgn(LT(M, pr).val) < 0)
            fmpz_sparse_vec_neg(&M->rows[pr], &M->rows[pr]);

        /* Reduce previous pivot rows incident to pivot column */
        for (i = nnz - npp; i < nnz; ++i)
            _fmpz_sparse_mat_with_transpose_gauss_elim(MT, pr, irows[i]);
        P[pr] = rank++;
    }

    /* Apply row permutation */
    fmpz_sparse_mat_permute_rows (M, P);

    flint_free(P);
    flint_free(irows);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(one);
    _fmpz_sparse_mat_with_transpose_clear(MT);
    return rank;
}
