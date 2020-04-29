/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

slong
fmpz_sparse_mat_howell_form_mod(fmpz_sparse_mat_t M, const fmpz_t mod)
{
    slong i, *P, rank = 0, remr = M->r;

    if (fmpz_sparse_mat_is_zero(M)) return 0;

    fmpz_sparse_mat_strong_echelon_form_mod(M, mod);
    P = flint_malloc(M->r*sizeof(*P));
    for (i = 0; i < M->r; ++i)
    {
        if (M->rows[i].nnz > 0) P[i] = rank++;
        else P[i] = --remr;
    }
    /* Apply row permutation */
    fmpz_sparse_mat_permute_rows (M, P);
    flint_free(P);
    return rank;
}

