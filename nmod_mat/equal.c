/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

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
            if (mat1->rows[i][j] != mat2->rows[i][j])
                return 0;
        }

        /*
        if (!_nmod_vec_equal(mat1->rows[i], mat2->rows[i], mat1->c))
            return 0;
        */
    }

    return 1;
}
