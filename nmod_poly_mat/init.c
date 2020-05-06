/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_init(nmod_poly_mat_t A, slong rows, slong cols, mp_limb_t n)
{
    if (rows > 0 && cols > 0)
    {
        slong i;

        A->entries = (nmod_poly_struct *) flint_malloc(flint_mul_sizes(rows, cols) * sizeof(nmod_poly_struct));
        A->rows = (nmod_poly_struct **) flint_malloc(rows * sizeof(nmod_poly_struct *));

        for (i = 0; i < rows * cols; i++)
            nmod_poly_init(A->entries + i, n);

        for (i = 0; i < rows; i++)
            A->rows[i] = A->entries + i * cols;

        A->r = rows;
        A->c = cols;
        A->modulus = n;
    }
    else
        memset(A, 0, sizeof(*A));
}
