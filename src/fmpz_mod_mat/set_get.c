/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod_mat.h"

/* Setters ********************************************************************/

void fmpz_mod_mat_set_entry(fmpz_mod_mat_t mat, slong i, slong j, const fmpz_t val)
{
    fmpz_set(fmpz_mat_entry(mat->mat, i, j), val);
}

void fmpz_mod_mat_set_fmpz_mat(fmpz_mod_mat_t A, const fmpz_mat_t B)
{
    fmpz_mat_set(A->mat, B);
    _fmpz_mod_mat_reduce(A);
}

void _fmpz_mod_mat_set_mod(fmpz_mod_mat_t mat, const fmpz_t n)
{
    fmpz_set(mat->mod, n);
}

void fmpz_mod_mat_set(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
{
    fmpz_set(B->mod, A->mod);
    fmpz_mat_set(B->mat, A->mat);
}

/* Getters ********************************************************************/

void fmpz_mod_mat_get_entry(fmpz_t x, const fmpz_mod_mat_t mat, slong i, slong j)
{
  fmpz_set(x, fmpz_mod_mat_entry(mat, i, j));
}

void fmpz_mod_mat_get_fmpz_mat(fmpz_mat_t A, const fmpz_mod_mat_t B)
{
    fmpz_mat_set(A, B->mat);
}
