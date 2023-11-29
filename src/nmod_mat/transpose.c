/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

void
nmod_mat_transpose(nmod_mat_t B, const nmod_mat_t A)
{
    mp_limb_t tmp;

    slong i, j;

    if (B->r != A->c || B->c != A->r)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_mat_transpose). Incompatible dimensions.\n");
    }

    if (A == B) /* In-place, guaranteed to be square */
    {
        for (i = 0; i < A->r - 1; i++)
            for (j = i + 1; j < A->c; j++)
            {
                tmp = A->rows[i][j];
                A->rows[i][j] = A->rows[j][i];
                A->rows[j][i] = tmp;
            }
    }
    else /* Not aliased; general case */
    {
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                B->rows[i][j] = A->rows[j][i];
    }
}
