/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

slong
nmod_sparse_mat_howell_form(nmod_sparse_mat_t M)
{
slong i, *P, rank = 0, remr = M->r;

    if (nmod_sparse_mat_is_zero(M)) return 0;

    nmod_sparse_mat_strong_echelon_form(M);
    P = flint_malloc(M->r*sizeof(*P));
    for (i = 0; i < M->r; ++i)
    {
        if (M->rows[i].nnz > 0) P[i] = rank++;
        else P[i] = --remr;
    }
    /* Apply row permutation */
    nmod_sparse_mat_permute_rows (M, P);
    flint_free(P);
    return rank;
}

