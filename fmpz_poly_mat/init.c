/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_init(fmpz_poly_mat_t A, slong rows, slong cols)
{
    slong i;

    if (rows != 0)
        A->rows = (fmpz_poly_struct **) flint_malloc(rows * sizeof(fmpz_poly_struct *));
    else
        A->rows = NULL;

    if (rows != 0 && cols != 0)
    {
        A->entries = (fmpz_poly_struct *) flint_malloc(flint_mul_sizes(rows, cols) * sizeof(fmpz_poly_struct));

        for (i = 0; i < rows * cols; i++)
            fmpz_poly_init(A->entries + i);

        for (i = 0; i < rows; i++)
            A->rows[i] = A->entries + i * cols;
    }
    else
    {
        A->entries = NULL;
        if (rows != 0)
        {
            for (i = 0; i < rows; i++)
                A->rows[i] = NULL;
        }
    }

    A->r = rows;
    A->c = cols;
}
