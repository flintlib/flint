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
#include "flint.h"
#include "fmpz_sparse_mat.h"

void fmpz_sparse_mat_from_entries(fmpz_sparse_mat_t M, slong * rows, slong * cols, fmpz * vals, slong nnz)
{
    slong r, i, j;
    for (r = i = 0; r < M->r; ++r, i = j)
    {
        for (j = i; j < nnz && rows[j]==r; ++j);
        fmpz_sparse_vec_from_entries(&M->rows[r], cols+i, vals+i, j-i);
    }
}
