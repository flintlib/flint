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

/* Elements must be ordered by row then col */
void nmod_sparse_mat_set_from_entries(nmod_sparse_mat_t mat, slong * rows, slong * cols, mp_limb_t * vals, slong nnz) 
{
    if(nnz==0) {
        nmod_sparse_mat_zero(mat);
        return;
    }
    slong i;
    mat->entries = flint_realloc(mat->entries, nnz * sizeof(*mat->entries));
    memset(mat->row_starts, 0, sizeof(*mat->row_starts));
    memset(mat->row_nnz, 0, sizeof(*mat->row_nnz));
    mat->c = 0;
    for(i = 0; i < nnz; i++) {
        if(i > 0  && rows[i] != rows[i-1]) {
            mat->row_starts[rows[i]] = i;
        }
        mat->entries[i].col = cols[i];
        mat->entries[i].val = vals[i];
        mat->row_nnz[rows[i]]++;
    }
    _nmod_sparse_mat_set_c(mat);
    mat->nnz = nnz;
}
