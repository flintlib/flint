/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

slong
fmpz_sparse_mat_hnf_modular(fmpz_sparse_mat_t M, const fmpz_t det)
{
    slong i, r, pr, pc, rank, remr, nnz;
    slong *P; /* Row permutation */
    slong *irows; /* Set of rows incident on a given column */
    slong nnp; /* Number of such rows which are not previous pivots */
    slong npp; /* Number of such rows which are previous pivots */
    fmpz_sparse_mat_with_transpose_t MT;
    hashmap_struct *hcol; /* Virtual L^t and one of its rows */
    fmpz_t g, a, b, one, rem_det;

    if (M->r == 0 || M->c == 0 || M->r != M->c || fmpz_is_zero(det)) return 0;
    fmpz_init(g);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init_set_ui(one, UWORD(1));
    fmpz_init_set(rem_det, det);

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

    irows = NULL;
    for (rank = pc = 0; pc < M->c; ++pc)
    {
        hcol = &MT->cols[pc]; nnz = hcol->num;
        if (!nnz) continue;
        irows = flint_realloc(irows, nnz*sizeof(*irows));
        pr = -1, nnp = 0, npp = 0;

        /* Find incident row pr which is not a previous pivot and has minimal leading term */
        /* Make pi_0 ... pi_{nnp-1} the other incident rows which are not previous pivots */
        /* and pi_{nnz-npp} ... pi_{nnz - 1} the incident rows which are previous pivots */
        for (i = 0; i < nnz; ++i)
        {
            r = hcol->keys[i];
            if (P[r] >= 0) 
                irows[hcol->num - (++npp)] = r;
            else if (pr >= 0 && fmpz_cmpabs(LT(M, r).val, LT(M, pr).val) >= 0) 
                irows[nnp++] = r;
            else
            {
                if (pr >= 0) irows[nnp++] = pr;
                pr = r;
            }
        }
        if (pr == -1) 
        {
            /* Set pivot col in some non-pivot row to rem_det */
            for (pr = 0; pr < M->r; ++pr)
                if (P[pr] == -1) break;
            fmpz_sparse_vec_set_entry(&M->rows[pr], pc, rem_det);
            fmpz_one(rem_det);
        }
        else
        {
            /* Eliminate non-pivot rows */
            for (i = 0; i < nnp; ++i)
                _fmpz_sparse_mat_with_transpose_gauss_elim_ext_mods(MT, pr, irows[i], rem_det);

            /* Minimize row modulo rem_det */
            fmpz_xgcd(g, a, b, LT(M, pr).val, rem_det);
        
            MT_FIX(MT, pr, 
            fmpz_sparse_vec_scalar_mul_fmpz(&M->rows[pr], &M->rows[pr], a);
            fmpz_sparse_vec_scalar_mods_fmpz(&M->rows[pr], &M->rows[pr], rem_det);
            if (M->rows[pr].nnz == 0 || LT(M, pr).ind != pc)
                fmpz_sparse_vec_set_entry(&M->rows[pr], pc, rem_det);
            );
            
            fmpz_divexact(rem_det, rem_det, g);
        }
        /* Reduce previous pivot rows */
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
    fmpz_clear(rem_det);
    _fmpz_sparse_mat_with_transpose_clear(MT);
    return M->r;
}
