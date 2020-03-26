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
void nmod_sparse_mat_append_row(nmod_sparse_mat_t mat, slong row, slong * cols, mp_limb_t * vals, slong nnz) 
{
    slong i;
    mat->entries = flint_realloc(mat->entries, (mat->nnz + nnz)* sizeof(*mat->entries));
    mat->row_starts[row] = mat->nnz;
    mat->row_nnz[row] = nnz;

    nmod_sparse_mat_entry_struct *e = mat->entries + mat->row_starts[row];
    for(i = 0; i < nnz; i++) {
        e[i].col = cols[i];
        e[i].val = vals[i];
        if(e[i].col >= mat->c) mat->c = e->col + 1;
    }
    mat->nnz += nnz;
}
