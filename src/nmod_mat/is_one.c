/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

int nmod_mat_is_one(const nmod_mat_t mat)
{
    slong i;

    if (mat->mod.n == 0 || mat->r == 0 || mat->c == 0)
        return 1;

    for (i = 0; i < mat->r; i++)
    {
        if (!_nmod_vec_is_zero(mat->rows[i], FLINT_MIN(mat->c, i)))
            return 0;

        if (i + 1 > mat->c)
            continue;

        if (mat->rows[i][i] != 1)
            return 0;

        if (!_nmod_vec_is_zero(mat->rows[i] + i + 1, mat->c - (i + 1)))
            return 0;
    }

    return 1;
}

