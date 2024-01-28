/*
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

int fmpz_mod_mat_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                                                        const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    slong i, rank, *perm;
    fmpz_mod_mat_t LU;
    int result;

    if (fmpz_mod_mat_is_empty(A, ctx))
        return 1;

    fmpz_mod_mat_init_set(LU, A, ctx);

    perm = flint_malloc(sizeof(slong) * A->r);
    for (i = 0; i < A->r; i++)
        perm[i] = i;

    rank = fmpz_mod_mat_lu(perm, LU, 1, ctx);

    if (rank == A->r)
    {
        fmpz_mod_mat_t PB;
        fmpz_mod_mat_window_init(PB, B, 0, 0, B->r, B->c, ctx);
        for (i = 0; i < A->r; i++)
            PB->rows[i] = B->rows[perm[i]];

        fmpz_mod_mat_solve_tril(X, LU, PB, 1, ctx);
        fmpz_mod_mat_solve_triu(X, LU, X, 0, ctx);

        fmpz_mod_mat_window_clear(PB, ctx);
        result = 1;
    }
    else
    {
        result = 0;
    }

    fmpz_mod_mat_clear(LU, ctx);
    flint_free(perm);

    return result;
}
