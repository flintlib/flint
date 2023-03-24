/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_mat.h"

int fmpz_mod_mat_equal(const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2)
{
    return fmpz_equal(mat1->mod, mat2->mod) && fmpz_mat_equal(mat1->mat, mat2->mat);
}

int fmpz_mod_mat_is_zero(const fmpz_mod_mat_t mat)
{
    return fmpz_mat_is_zero(mat->mat);
}

int fmpz_mod_mat_is_one(const fmpz_mod_mat_t mat)
{
    return fmpz_is_one(mat->mod) || fmpz_mat_is_one(mat->mat);
}
