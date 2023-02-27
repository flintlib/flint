/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_mat.h"

void
d_mat_randtest(d_mat_t mat, flint_rand_t state, slong minexp, slong maxexp)
{
    slong r, c, i, j;

    r = mat->r;
    c = mat->c;

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            d_mat_entry(mat, i, j) = d_randtest_signed(state, minexp, maxexp);
}
