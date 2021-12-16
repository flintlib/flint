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

int fmpz_mod_mat_inv(fmpz_mod_mat_t B, fmpz_mod_mat_t A)
{
    fmpz_mod_mat_t I;
    slong i, dim;
    int result;

    dim = A->mat->r;

    switch (dim)
    {
        case 0:
            result = 1;
            break;

        case 1:
            if (fmpz_is_zero(fmpz_mod_mat_entry(A, 0, 0)))
            {
                result = 0;
            }
            else
            {
                fmpz_mod_ctx_t ctx;
                fmpz_mod_ctx_init(ctx, A->mod);
                fmpz_mod_inv(fmpz_mod_mat_entry(B, 0, 0),
                             fmpz_mod_mat_entry(A, 0, 0), ctx);
                fmpz_mod_ctx_clear(ctx);
                result = 1;
            }
            break;

        default:
            fmpz_mod_mat_init(I, dim, dim, A->mod);
            for (i = 0; i < dim; i++)
                fmpz_one(fmpz_mod_mat_entry(I, i, i));
            result = fmpz_mod_mat_solve(B, A, I);
            fmpz_mod_mat_clear(I);
    }

    return result;
}

