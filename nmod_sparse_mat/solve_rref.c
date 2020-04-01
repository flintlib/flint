/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "nmod_sparse_mat.h"

int nmod_sparse_mat_solve_rref(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b)
{
    int good = 1;
    slong i;
    nmod_sparse_mat_t AB;
    nmod_sparse_vec_struct *row;
    nmod_sparse_entry_struct *le, *re;
    nmod_sparse_mat_init(AB, A->r, A->c, A->mod);
    nmod_sparse_mat_set(AB, A);
    nmod_sparse_mat_append_col(AB, b);
    AB->c = A->c;
    nmod_sparse_mat_rref(AB);
    AB->c = A->c+1;

    memset(x, 0, A->c*sizeof(*x));
    for (i = 0; i < A->r; ++i) 
    {
        row = &AB->rows[i];
        if (row->nnz == 0) continue;
        le = &row->entries[0];
        re = &row->entries[row->nnz-1];
        /* If any row has leading col A->c, system is not solvable */
        if (le->ind==A->c) {good = 0; break;}
        /* Otherwise, x[lc] = lagging value */
        if (re->ind==A->c) {x[le->ind] = re->val;}
    }
    nmod_sparse_mat_clear(AB);
}
