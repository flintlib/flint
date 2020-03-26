/*
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

/*Function to convert a square matrix to an identity matrix
    The matrix is assumed to be a square one*/
void
nmod_sparse_mat_one(nmod_sparse_mat_t mat)
{
    mat->c = mat->nnz = mat->r;
    mat->entries = flint_realloc(mat->entries, mat->r*sizeof(*mat->entries));
    slong i;
    for(i = 0; i < mat->r; i++) {
        mat->entries[i].val = 1;
        mat->entries[i].col = i;
        mat->row_starts[i] = i;
        mat->row_nnz[i] = 1;
    }
}
