/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_mat.h"

/* Setters ********************************************************************/

void fmpz_mod_mat_set_entry(fmpz_mod_mat_t mat, slong i, slong j, const fmpz_t val, const fmpz_mod_ctx_t ctx)
{
    fmpz_set(fmpz_mat_entry(mat, i, j), val);
}

void fmpz_mod_mat_set(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_set(B, A);
}

/* Getters ********************************************************************/

void fmpz_mod_mat_get_entry(fmpz_t x, const fmpz_mod_mat_t mat, slong i, slong j, const fmpz_mod_ctx_t ctx)
{
    fmpz_set(x, fmpz_mod_mat_entry(mat, i, j));
}

void fmpz_mod_mat_get_fmpz_mat(fmpz_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_set(A, B);
}
