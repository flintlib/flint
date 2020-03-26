/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "nmod_sparse_mat.h"

void
nmod_sparse_mat_concat_horizontal(nmod_sparse_mat_t res, const nmod_sparse_mat_t mat1, const nmod_sparse_mat_t mat2)
{
    if(res->r == 0) return;
    if(mat1->c == 0) {
        nmod_sparse_mat_set(res, mat2);
        return;
    }
    if(mat2->c == 0) {
        nmod_sparse_mat_set(res, mat1);
        return;
    }
    slong i, j;
    res->c = mat1->c + mat2->c;
    res->nnz = mat1->nnz + mat2->nnz;
    res->entries = flint_realloc(res->entries, res->nnz * sizeof(*res->entries));
    
    for (i = 0; i < res->r; i++)
    {
        res->row_starts[i] = mat1->row_starts[i] + mat2->row_starts[i];
        nmod_sparse_mat_entry_struct *e = res->entries + res->row_starts[i];
        nmod_sparse_mat_entry_struct *e1 = mat1->entries + mat1->row_starts[i];
        nmod_sparse_mat_entry_struct *e2 = mat2->entries + mat2->row_starts[i];
        for(j = 0; j < mat1->row_nnz[i]; ++j, ++e, ++e1) {
            e->col = e1->col;
            e->val = e1->val;
        }
        for(j = 0; j < mat2->row_nnz[i]; ++j, ++e, ++e2) {
            e->col = e2->col + mat1->c;
            e->val = e2->val;
        }
        res->row_nnz[i] = mat1->row_nnz[i] + mat2->row_nnz[i];
    }
}
