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

/* Assumes B->r == A->r */
void
nmod_sparse_mat_set(nmod_sparse_mat_t B, const nmod_sparse_mat_t A)
{
    if (B == A)
        return;
    B->c = A->c;
    B->nnz = A->nnz;
    if(B->nnz == 0) {
        flint_free(B->entries);
        B->entries = NULL;
    } else {
        B->entries = flint_realloc(B->entries, A->nnz * sizeof(*B->entries));
        memcpy(B->entries, A->entries, A->nnz * sizeof(*B->entries));
    }
    memcpy(B->row_starts, A->row_starts, A->r * sizeof(*B->row_starts));
    memcpy(B->row_nnz, A->row_nnz, A->r * sizeof(*B->row_nnz));
}
