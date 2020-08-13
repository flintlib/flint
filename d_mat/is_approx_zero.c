/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_mat.h"

int
d_mat_is_approx_zero(const d_mat_t mat, double eps)
{
    slong j;

    if (mat->r == 0 || mat->c == 0)
        return 1;

    for (j = 0; j < mat->r; j++)
    {
        if (!_d_vec_is_approx_zero(mat->rows[j], mat->c, eps))
            return 0;
    }

    return 1;
}
