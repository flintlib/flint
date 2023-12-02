/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, slong rank,
                  flint_bitcnt_t bits)
{
    slong i;
    fmpz * diag;

    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_randrank). Impossible rank.\n");
    }

    diag = _fmpz_vec_init(rank);
    for (i = 0; i < rank; i++)
        fmpz_randtest_not_zero(&diag[i], state, bits);

    fmpz_mat_randpermdiag(mat, state, diag, rank);

    _fmpz_vec_clear(diag, rank);
}
