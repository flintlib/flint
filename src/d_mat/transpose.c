/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_mat.h"

void
d_mat_transpose(d_mat_t B, const d_mat_t A)
{
    slong ii, jj, i, j, blocksize;
    blocksize = 64 / sizeof(double);

    if (B->r != A->c || B->c != A->r)
    {
        flint_throw(FLINT_ERROR, "Exception (d_mat_transpose). Incompatible dimensions.\n");
    }

    if (B == A)
    {
        d_mat_t t;
        d_mat_init(t, A->r, A->c);
        d_mat_transpose(t, A);
        d_mat_swap_entrywise(B, t);
        d_mat_clear(t);
        return;
    }

    for (ii = 0; ii < B->r; ii += blocksize)
        for (jj = 0; jj < B->c; jj += blocksize)
            for (i = ii; i < FLINT_MIN(ii + blocksize, B->r); i++)
                for (j = jj; j < FLINT_MIN(jj + blocksize, B->c); j++)
                    d_mat_entry(B, i, j) = d_mat_entry(A, j, i);

}
