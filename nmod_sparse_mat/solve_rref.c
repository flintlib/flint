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

int nmod_sparse_mat_solve_rref(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b)
{
    int good = 1;
    slong i;
    nmod_sparse_mat_t Mb;
    nmod_sparse_vec_struct *row;
    nmod_sparse_entry_struct *le, *re;
    if (_nmod_vec_is_zero(b, M->c))
    {
        _nmod_vec_zero(x, M->c);
        return 1;
    }

    nmod_sparse_mat_init(Mb, M->r, M->c, M->mod);
    nmod_sparse_mat_set(Mb, M);
    nmod_sparse_mat_append_col(Mb, b);
    Mb->c = M->c;
    nmod_sparse_mat_rref(Mb);
    Mb->c = M->c+1;

    memset(x, 0, M->c*sizeof(*x));
    for (i = 0; i < M->r; ++i) 
    {
        row = &Mb->rows[i];
        if (row->nnz == 0) continue;
        le = &row->entries[0];
        re = &row->entries[row->nnz-1];
        /* If any row has leading col M->c, system is not solvable */
        if (le->ind==M->c) {good = 0; break;}
        /* Otherwise, x[lc] = lagging value */
        if (re->ind==M->c) {x[le->ind] = re->val;}
    }
    nmod_sparse_mat_clear(Mb);
    return good;
}
