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
fmpz_mat_neg(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong i;

    if (res->c < 1)
        return;

    for (i = 0; i < res->r; i++)
        _fmpz_vec_neg(fmpz_mat_row(res, i), fmpz_mat_row(mat, i), res->c);
}
