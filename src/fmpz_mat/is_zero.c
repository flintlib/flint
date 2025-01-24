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

int
fmpz_mat_is_zero(const fmpz_mat_t mat)
{
    slong j;

    if (mat->r == 0 || mat->c == 0)
        return 1;

    for (j = 0; j < mat->r; j++)
    {
        if (!_fmpz_vec_is_zero(fmpz_mat_row(mat, j), mat->c))
            return 0;
    }

    return 1;
}

int
fmpz_mat_is_zero_row(const fmpz_mat_t mat, slong i)
{
    return _fmpz_vec_is_zero(fmpz_mat_row(mat, i), mat->c);
}
