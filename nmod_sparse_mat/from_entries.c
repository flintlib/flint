/*
<<<<<<< HEAD
<<<<<<< HEAD
    Copyright (C) 2011 Fredrik Johansson
=======
    Copyright (C) 2010 Fredrik Johansson
>>>>>>> Initial code for sparse matrices mod limb size integers, just construction and arithmetic for starters
=======
    Copyright (C) 2011 Fredrik Johansson
>>>>>>> Added sparse vector class to nmod, changed sparse matrix class to use it for underlying, added (untested) LU decomposition

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

void nmod_sparse_mat_from_entries(nmod_sparse_mat_t M, slong * rows, slong * cols, mp_limb_t * vals, slong nnz)
{
    slong r, i, j;
    for (r = i = 0; r < M->r; ++r, i = j)
    {
        M->rows[r].nnz = 0;
        for (j = i; j < nnz && rows[j]==r; ++j);
        nmod_sparse_vec_from_entries(&M->rows[r], cols+i, vals+i, j-i);
    }
}
