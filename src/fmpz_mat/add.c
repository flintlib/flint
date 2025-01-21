/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_add(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2)
{
    slong i;

    if (res->c < 1)
        return;

    for (i = 0; i < res->r; i++)
        _fmpz_vec_add(fmpz_mat_row(res, i),
            fmpz_mat_row(mat1, i), fmpz_mat_row(mat2, i), res->c);
}
