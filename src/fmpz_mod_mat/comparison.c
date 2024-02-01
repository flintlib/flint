/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_mat.h"

int fmpz_mod_mat_equal(const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)
{
    return fmpz_mat_equal(mat1, mat2);
}

int fmpz_mod_mat_is_zero(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    return fmpz_mat_is_zero(mat);
}

int fmpz_mod_mat_is_one(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    return fmpz_is_one(ctx->n) || fmpz_mat_is_one(mat);
}
