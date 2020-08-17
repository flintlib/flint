/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_zero(fmpz_mat_t mat)
{
    slong i;

    if (mat->c < 1)
        return;

    for (i = 0; i < mat->r; i++)
        _fmpz_vec_zero(mat->rows[i], mat->c);
}
