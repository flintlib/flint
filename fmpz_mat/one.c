/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_one(fmpz_mat_t mat)
{
    slong i, n;

    fmpz_mat_zero(mat);
    n = FLINT_MIN(mat->r, mat->c);

    for (i = 0; i < n; i++)
        fmpz_one(fmpz_mat_entry(mat, i, i));
}
