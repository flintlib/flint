/*
    Copyright (C) 2015 Tommy Hofmann 

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"




/*  Multiply by unit to minimize leading term (to a factor of N) */
static void scale_row(fmpz_sparse_mat_with_transpose_t MT, slong r, const fmpz_t mod)
{
    fmpz_t g, a, b, n;
    if (!fmpz_is_one(LT(MT->M, r).val))
    {
        fmpz_init(g);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_xgcd(g, a, b, LT(MT->M, r).val, mod);
        if (!fmpz_is_one(g)) /* Need to lift a = (M[pr][pc]/g)^-1 mod N/g to a unit modulo N */
        {
            fmpz_divexact(b, mod, g);
            fmpz_init_set(n, mod);
            while (!fmpz_is_one(g))
            {
                fmpz_gcd(g, a, n);
                fmpz_divexact(n, n, g);
            }
            fmpz_addmul(a, n, b); 
            fmpz_mod(a, a, mod);
            fmpz_clear(n);
        }
        MT_FIX(MT, r, 
        fmpz_sparse_vec_scalar_mul_fmpz(&MT->M->rows[r], &MT->M->rows[r], a);
        fmpz_sparse_vec_scalar_mod_fmpz(&MT->M->rows[r], &MT->M->rows[r], mod);
        );
        fmpz_clear(g);
        fmpz_clear(a);
        fmpz_clear(b);
    }
}

slong
fmpz_sparse_mat_strong_echelon_form_mod(fmpz_sparse_mat_t M, const fmpz_t mod)
{
    slong i, r, c, pr, pc, rank, remr, nzrows;
    slong *P; /* Row permutation */
    slong *Pr; /* Map from column to associated pivot row */
    slong *irows; /* Set of rows incident on a given column */
    slong *zrows; /* Set of empty rows */
    slong nnp; /* Number of such rows which are not pivots */
    slong npp; /* Number of such rows which are not pivots */
    fmpz_sparse_mat_with_transpose_t MT;
    hashmap_struct *hcol; /* Virtual L^t and one of its rows */
    fmpz_t q, zero;
    fmpz_sparse_vec_t zero_vec;

    if (fmpz_sparse_mat_is_zero(M)) return 0;
    fmpz_init(zero);
    fmpz_init(q);

    /* Final object must have at least as many rows as columns */
    fmpz_sparse_vec_init(zero_vec);
    while (M->r < M->c) fmpz_sparse_mat_append_row(M, zero_vec);

    /* Need extra row to deal with final elimination */
    fmpz_sparse_mat_append_row(M, zero_vec);

    /* Initialize data structure to hold copy of incident rows of a given column */
    irows = flint_malloc(M->r*sizeof(*irows));
    zrows = flint_malloc(M->r*sizeof(*zrows));

    fmpz_sparse_mat_scalar_mod_fmpz(M, M, mod);

    /* Set up permutation */
    P = flint_malloc(M->r*sizeof(*P));
    remr = M->r, nzrows = 0;
    for (r = 0; r < M->r; ++r) 
    {
        if (!M->rows[r].nnz) P[r] = --remr, zrows[nzrows++] = r; 
        else P[r] = -1;
    }

    /* Set up pivot row mapping */
    Pr = flint_malloc(M->c*sizeof(*Pr));
    for (c = 0; c < M->c; ++c)
        Pr[c] = -1;

    /* Construct virtual transpose */
    _fmpz_sparse_mat_with_transpose_init(MT, M);
    
    for (pc = 0; pc < M->c; ++pc)
    {
        hcol = &MT->cols[pc];
        if (!hcol->num) continue;
        pr = -1, nnp = 0;

        /* Find incident row pr which is not a previous pivot and has minimal leading term */
        /* Make pi_0 ... pi_{nnp-1} the other incident rows which are not previous pivots */
        for (i = 0; i < hcol->num; ++i)
        {
            r = hcol->keys[i];
            if (P[r] >= 0) continue;
            else if (pr >= 0 && fmpz_cmpabs(LT(M, r).val, LT(M, pr).val) >= 0) 
                irows[nnp++] = r;
            else
            {
                if (pr >= 0) irows[nnp++] = pr;
                pr = r;
            }
        }
        if (pr == -1) continue; /* Cannot perform elimination on this column (yet) */
        scale_row(MT, pr, mod);
        Pr[pc] = pr, P[pr] = pc;

        /* Eliminate non-pivot rows */
        for (i = 0; i < nnp; ++i)
        {
            r = irows[i];
            _fmpz_sparse_mat_with_transpose_gauss_elim_ext_mod(MT, pr, r, mod);
            if (M->rows[r].nnz == 0) P[r] = --remr, zrows[nzrows++] = r;
        }
    }

    /* Reduce upwards, and deal with non-unit leading terms */
    for (pc = 0; pc < M->c; pc++)
    {
        if (Pr[pc] == -1) {P[zrows[--nzrows]] = pc; continue;} /* No pivot for this column */
        hcol = &MT->cols[pc];
        pr = Pr[pc];

        /* Reduce previous pivot rows */
        for (i = npp = 0; i < hcol->num; ++i)
            if (hcol->keys[i] != pr)
                irows[npp++] = hcol->keys[i];
        for (i = 0; i < npp; ++i)
        {
            _fmpz_sparse_mat_with_transpose_gauss_elim_mod(MT, pr, irows[i], mod);
        }
        if (fmpz_is_one(LT(M, pr).val)) continue;

        /* Obtain (possibly) new basis element by scalar multiplicition */
        r = zrows[--nzrows];
        fmpz_divexact(q, mod, LT(M, pr).val);
        MT_FIX(MT, r,
        fmpz_sparse_vec_scalar_mul_fmpz(&M->rows[r], &M->rows[pr], q);
        fmpz_sparse_vec_scalar_mod_fmpz(&M->rows[r], &M->rows[r], mod);
        );

        while (!fmpz_sparse_vec_is_zero(&M->rows[r]))
        {
            c = M->rows[r].entries[0].ind;
            /* If no previous pivot row exists, use this row */
            if (Pr[c] == -1) 
            {
                scale_row(MT, r, mod);
                Pr[c] = r; P[r] = c; 
                break;
            }

            /* Otherwise, eliminate c using existing pivot */
            _fmpz_sparse_mat_with_transpose_gauss_elim_ext_mod(MT, Pr[c], r, mod);
            if (!fmpz_sparse_vec_is_zero(&M->rows[r]) && LT(M, r).ind == c) flint_abort();
        }
        /* If row fully eliminated, add back to empty stock */

        if (fmpz_sparse_vec_is_zero(&M->rows[r])) nzrows++;
    
    }
    fmpz_sparse_mat_permute_rows(M, P);
    fmpz_clear(zero);
    fmpz_clear(q);
    flint_free(irows);
    flint_free(zrows);
    fmpz_sparse_vec_clear(zero_vec);
    M->r -= 1;
    M->rows = realloc(M->rows, M->r*sizeof(*M->rows));
    _fmpz_sparse_mat_with_transpose_clear(MT);
    flint_free(P);
    flint_free(Pr);
    return rank;
}

