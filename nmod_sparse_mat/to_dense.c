/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

void nmod_sparse_mat_to_dense(nmod_mat_t mat, const nmod_sparse_mat_t src)
{
    if(mat->r == 0 || mat->c == 0) return;
    slong i, j, k=0;
    memset(mat->entries, 0, mat->r * mat->c * sizeof(*mat->entries));
    for(i=0; i<src->r; ++i) {
        for(j=0; j<src->row_nnz[i]; ++j, ++k) {
            mat->rows[i][src->entries[k].col] = src->entries[k].val;
        }
    }
}
