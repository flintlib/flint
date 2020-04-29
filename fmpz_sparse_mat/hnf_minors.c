/*
    Copyright (C) 2014 Alex J. Best
    Copyright (C) 2017 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

/*
  This is the algorithm of Kannan, Bachem, "Polynomial algorithms for computing
  the Smith and Hermite normal forms of an integer matrix", Siam J. Comput.,
  Vol. 8, No. 4, pp. 499-507.
*/
slong
fmpz_sparse_mat_hnf_minors(fmpz_sparse_mat_t M)
{

    slong i, r, c, pr, pc, rank, remr, nnz;
    slong *P; /* Row permutation */
    slong *Pr; /* Assignment of pivot row to each column */
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
    
    /* Set up permutation and its "inverse" */
    P = flint_malloc(M->r*sizeof(*P));
    Pr = flint_malloc(M->c*sizeof(*Pr));
    remr = M->r;
    for (r = 0; r < M->r; ++r) 
    {
        if (!M->rows[r].nnz) P[r] = --remr;
        else P[r] = -1;
    }

    for (rank = pc = 0; pc < M->c; ++pc)
    {
        Pr[pc] = -1;

        /* Find first non-empty row which is not a previous pivot and eliminate cols up to pc */
        for (r = 0; r < M->r; ++r)
        {
            if (P[r] >= 0) continue;

            /* Reduce r by previous pivot rows */
            while ((c = M->rows[r].entries[0].ind) < pc)
            {
                _fmpz_sparse_mat_with_transpose_gauss_elim_ext(MT, Pr[c], r);
                if (M->rows[r].nnz == 0) {P[r] = --remr; break;}
            }
            if (M->rows[r].nnz > 0 && M->rows[r].entries[0].ind == pc) break;
        }
        if (r == M->r) continue; /* No viable pivot for column */
        pr = r;

        if (fmpz_sgn(LT(M, r).val) < 0)
            fmpz_sparse_vec_neg(&M->rows[r], &M->rows[r]);

        /* Use this row to reduce previous pivot rows */
        hcol = &MT->cols[pc]; nnz = hcol->num;
        for (i = 0; i < nnz; ++i)
        {
            r = hcol->keys[i];
            if (r != pr && P[r] >= 0) 
            _fmpz_sparse_mat_with_transpose_gauss_elim(MT, pr, r);
        }
        Pr[pc] = pr;
        P[pr] = rank++;
    }

    /* Deal with any remaining rows */
    for (r = 0; r < M->r; ++r)
    {
        if (P[r] >= 0) continue;

        /* Reduce r by previous pivot rows */
        while ((c = M->rows[r].entries[0].ind) < M->c)
        {
            _fmpz_sparse_mat_with_transpose_gauss_elim_ext(MT, Pr[c], r);
            if (M->rows[r].nnz == 0) {P[r] = --remr; break;}
        }
    }

    /* Since pivot rows were modified, need to re-reduce previous pivot rows */
    for (pc = 0; pc < M->c; ++pc)
    {
        if (Pr[pc] == -1) continue; /* No pivot for this column */
        pr = Pr[pc];
        hcol = &MT->cols[pc]; nnz = hcol->num;
        for (i = 0; i < nnz; ++i)
        {
            r = hcol->keys[i];
            if (r == pr) continue;

            /* All other incident rows must be previous pivots */
            _fmpz_sparse_mat_with_transpose_gauss_elim(MT, pr, r);
        }
    }

    /* Apply row permutation */
    fmpz_sparse_mat_permute_rows (M, P);

    flint_free(P);
    flint_free(Pr);
    _fmpz_sparse_mat_with_transpose_clear(MT);
    fmpz_clear(g);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(one);
    return rank;
}
