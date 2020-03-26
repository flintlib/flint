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
nmod_sparse_mat_concat_vertical(nmod_sparse_mat_t res, const nmod_sparse_mat_t mat1, const nmod_sparse_mat_t mat2)
{
    slong i;
    
    res->c = (mat1->c > mat2->c)?(mat1->c):(mat2->c);
    res->nnz = mat1->nnz + mat2->nnz;
    res->entries = flint_realloc(res->entries, res->nnz * sizeof(*res->entries));
    memcpy(res->entries, mat1->entries, mat1->nnz * sizeof(*res->entries));
    memcpy(res->entries + mat1->nnz, mat2->entries, mat2->nnz * sizeof(*res->entries));
    memcpy(res->row_starts, mat1->row_starts, mat1->r * sizeof(*res->row_starts));
    memcpy(res->row_starts+mat1->r, mat2->row_starts, mat2->r * sizeof(*res->row_starts));
    memcpy(res->row_nnz, mat1->row_nnz, mat1->r*sizeof(*res->row_nnz));
    memcpy(res->row_nnz+mat1->r, mat2->row_nnz, mat2->r*sizeof(*res->row_nnz));
    
    for (i = 0; i < mat2->r; i++)
    	res->row_starts[mat1->r + i] += mat1->nnz;
}
