/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "templates.h"

int TEMPLATE(T, sparse_mat_solve_rref) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, const TEMPLATE(T, ctx_t) ctx)
{
    int good = 1;
    slong i;
    TEMPLATE(T, sparse_mat_t) Mb;
    TEMPLATE(T, sparse_vec_struct) *row;
    TEMPLATE(T, sparse_entry_struct) *le, *re;
    if (_TEMPLATE(T, vec_is_zero) (b, M->r, ctx))
    {
        _TEMPLATE(T, vec_zero) (x, M->c, ctx);
        return 1;
    }

    TEMPLATE(T, sparse_mat_init) (Mb, M->r, M->c, ctx);
    
    TEMPLATE(T, sparse_mat_set) (Mb, M, ctx);
    TEMPLATE(T, sparse_mat_append_col) (Mb, b, ctx);
    Mb->c = M->c;
    TEMPLATE(T, sparse_mat_rref) (Mb, ctx);
    Mb->c = M->c+1;

    _TEMPLATE(T, vec_zero) (x, M->c, ctx);
    for (i = 0; i < M->r; ++i) 
    {
        row = &Mb->rows[i];
        if (row->nnz == 0) continue;
        le = &row->entries[0];
        re = &row->entries[row->nnz-1];
        /* If any row has leading col M->c, system is not solvable */
        if (le->ind==M->c) {good = 0; break;}
        /* Otherwise, x[lc] = lagging value */
        if (re->ind==M->c) TEMPLATE(T, set) (&x[le->ind], re->val, ctx);
    }
    TEMPLATE(T, sparse_mat_clear) (Mb, ctx);
    return good;
}

#endif
