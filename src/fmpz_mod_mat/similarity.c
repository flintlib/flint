/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void fmpz_mod_mat_similarity(fmpz_mod_mat_t A, slong r, fmpz_t d)
{
    slong n = A->mat->r, i, j;
    fmpz_t t;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init(ctx, A->mod);

    fmpz_init(t);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < r - 1; j++)
        {
            fmpz_mod_mul(t, fmpz_mod_mat_entry(A, i, r), d, ctx);
            fmpz_mod_add(fmpz_mod_mat_entry(A, i, j), fmpz_mod_mat_entry(A, i, j), t, ctx);
        }

        for (j = r + 1; j < n; j++)
        {
            fmpz_mod_mul(t, fmpz_mod_mat_entry(A, i, r), d, ctx);
            fmpz_mod_add(fmpz_mod_mat_entry(A, i, j), fmpz_mod_mat_entry(A, i, j), t, ctx);
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < r - 1; j++)
        {
            fmpz_mod_mul(t, fmpz_mod_mat_entry(A, j, i), d, ctx);
            fmpz_mod_sub(fmpz_mod_mat_entry(A, r, i), fmpz_mod_mat_entry(A, r, i), t, ctx);
        }

        for (j = r + 1; j < n; j++)
        {
            fmpz_mod_mul(t, fmpz_mod_mat_entry(A, j, i), d, ctx);
            fmpz_mod_sub(fmpz_mod_mat_entry(A, r, i), fmpz_mod_mat_entry(A, r, i), t, ctx);
        }
    }

    fmpz_clear(t);

    fmpz_mod_ctx_clear(ctx);
}

