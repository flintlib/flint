/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

int
nmod_mat_equal(const nmod_mat_t mat1, const nmod_mat_t mat2)
{
    slong i, j;

    if (mat1->r != mat2->r || mat1->c != mat2->c)
        return 0;

    if (mat1->r == 0 || mat1->c == 0)
        return 1;

    for (i = 0; i < mat1->r; i++)
    {
        for (j = 0; j < mat1->c; j++)
        {
            if (!_nmod_vec_equal(nmod_mat_entry_ptr(mat1, i, 0), nmod_mat_entry_ptr(mat2, i, 0), mat1->c))
                return 0;
        }
    }

    return 1;
}
