/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "hashmap.h"
#include "fmpz_sparse_mat.h"

slong
fmpz_sparse_mat_hnf_classical(fmpz_sparse_mat_t M)
{
    slong i, r, pr, pc, rank, remr, next_pr, nnz;
    slong *P; /* Row permutation */
    slong *irows; /* Set of row incident on a given column */
    slong nnp, next_nnp; /* Number of such rows which are not previous pivots */
    slong npp; /* Number of such rows which are previous pivots */
    fmpz_sparse_mat_with_transpose_t MT;
    hashmap_struct *hcol;

    if (M->r == 0 || M->c == 0) return 0;
    
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
        pr = next_pr = -1, next_nnp = nnz, npp = 0;
        
        do
        {
            pr = next_pr, next_pr = -1;
            nnp = next_nnp, next_nnp = 0;
            for (i = 0; i < nnp; ++i)
            {
                r = (pr == -1) ? hcol->keys[i] : irows[i];
                if (pr >= 0)  /* Reduce row r by row pr */
                {
                    if (_fmpz_sparse_mat_with_transpose_gauss_elim(MT, pr, r))
                    {
                        if (M->rows[r].nnz == 0) P[r] = --remr;
                        continue;
                    }
                }
                else if (P[r] >= 0) /* Record incident previous pivot rows at end of irows */ 
                {
                    irows[hcol->num - (++npp)] = r; 
                    continue;
                }

                /* Either add r (back) to irows or make it the next pivot row */

                if (next_pr >= 0 && fmpz_cmpabs(LT(M, r).val, LT(M, next_pr).val) >= 0) 
                {
                    irows[next_nnp++] = r; 
                    continue;
                }
                if (next_pr >= 0) irows[next_nnp++] = next_pr;
                next_pr = r;
            }
            if (pr != -1) irows[next_nnp++] = pr;
        } while (next_pr != -1); /* Stop when no more reduction to be done */
        if (pr == -1) continue; /* No incident rows which are not previous pivots */

        if (fmpz_sgn(LT(M, pr).val) < 0)
            fmpz_sparse_vec_neg(&M->rows[pr], &M->rows[pr]);

        /* Now use row pr to reduce the previous pivot rows */
        for (i = nnz - npp; i < nnz; ++i)
            _fmpz_sparse_mat_with_transpose_gauss_elim(MT, pr, irows[i]);
        P[pr] = rank++;
    }

    /* Apply row permutation */
    fmpz_sparse_mat_permute_rows (M, P);

    flint_free(P);
    flint_free(irows);
    _fmpz_sparse_mat_with_transpose_clear(MT);
    return rank;
}
