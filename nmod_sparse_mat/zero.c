/*
    Copyright (C) 2011 Fredrik Johansson

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

void
nmod_sparse_mat_zero(nmod_sparse_mat_t mat)
{
    memset(mat->row_nnz, 0, mat->r*sizeof(*mat->row_nnz));
    memset(mat->row_starts, 0, mat->r*sizeof(*mat->row_starts));
    flint_free(mat->entries);
    mat->entries = NULL;
    mat->nnz = 0;
}
