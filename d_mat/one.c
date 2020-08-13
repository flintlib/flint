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

void
d_mat_one(d_mat_t mat)
{
    slong i, n;

    d_mat_zero(mat);
    n = FLINT_MIN(mat->r, mat->c);

    for (i = 0; i < n; i++)
        d_mat_entry(mat, i, i) = 1;
}
