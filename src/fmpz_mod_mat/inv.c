/*
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

int fmpz_mod_mat_inv(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_mat_t I;
    slong dim;
    int result;

    dim = A->r;

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
                fmpz_mod_inv(fmpz_mod_mat_entry(B, 0, 0),
                             fmpz_mod_mat_entry(A, 0, 0), ctx);
                result = 1;
            }
            break;

        default:
            fmpz_mod_mat_init(I, dim, dim, ctx);
            fmpz_mod_mat_one(I, ctx);
            result = fmpz_mod_mat_solve(B, A, I, ctx);
            fmpz_mod_mat_clear(I, ctx);
    }

    return result;
}
