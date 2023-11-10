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

void fmpz_mod_mat_swap_rows(fmpz_mod_mat_t mat, slong * perm, slong r, slong s)
{
    fmpz_mat_swap_rows(mat->mat, perm, r, s);
}

void fmpz_mod_mat_invert_rows(fmpz_mod_mat_t mat, slong * perm)
{
    fmpz_mat_invert_rows(mat->mat, perm);
}

void fmpz_mod_mat_swap_cols(fmpz_mod_mat_t mat, slong * perm, slong r, slong s)
{
    fmpz_mat_swap_cols(mat->mat, perm, r, s);
}

void fmpz_mod_mat_invert_cols(fmpz_mod_mat_t mat, slong * perm)
{
    fmpz_mat_invert_cols(mat->mat, perm);
}
