/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson

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

void
nmod_sparse_mat_clear(nmod_sparse_mat_t mat)
{
    if(mat->entries)
    {
        flint_free(mat->entries);
    }
    if(mat->row_starts) {
        flint_free(mat->row_starts);
    }
    if(mat->row_nnz) {
        flint_free(mat->row_nnz);
    }
}
