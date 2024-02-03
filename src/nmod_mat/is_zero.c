/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

int
nmod_mat_is_zero(const nmod_mat_t mat)
{
    slong j;

    if (mat->r == 0 || mat->c == 0)
        return 1;

    for (j = 0; j < mat->r; j++)
    {
        if (!_nmod_vec_is_zero(mat->rows[j], mat->c))
            return 0;
    }

    return 1;
}
