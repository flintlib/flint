/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_randrank(fmpz_mod_mat_t mat, flint_rand_t state, slong rank, const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz * diag;

    if (rank < 0 || rank > fmpz_mod_mat_nrows(mat, ctx) || rank > fmpz_mod_mat_ncols(mat, ctx))
        flint_throw(FLINT_ERROR, "Impossible rank in %s\n", __func__);

    diag = _fmpz_vec_init(rank);

    for (i = 0; i < rank; i++)
    {
        fmpz_randm(diag + i, state, ctx->n);
        if (fmpz_is_zero(diag + i))
            fmpz_one(diag + i);
    }

    fmpz_mat_randpermdiag(mat, state, diag, rank);

    _fmpz_vec_clear(diag, rank);
}
