/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010,2011 Fredrik Johansson

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
nmod_sparse_mat_init(nmod_sparse_mat_t mat, slong rows, mp_limb_t n)
{
    memset(mat, 0, sizeof(*mat));
    mat->r = rows;
    mat->c = 0;
    if(rows > 0) {
        mat->row_starts = flint_calloc(rows, rows*sizeof(*mat->row_starts));
        mat->row_nnz = flint_calloc(rows, rows*sizeof(*mat->row_nnz));
    }
    _nmod_sparse_mat_set_mod(mat, n);
}
