/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz_mod_mat.h"

void fmpz_mod_mat_concat_horizontal(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2)
{
    fmpz_mat_concat_horizontal(res->mat, mat1->mat, mat2->mat);
}

void fmpz_mod_mat_concat_vertical(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2)
{
    fmpz_mat_concat_vertical(res->mat, mat1->mat, mat2->mat);
}
