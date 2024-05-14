/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "long_extras.h"
#include "fmpq_mat.h"

void fmpq_mat_init(fmpq_mat_t mat, slong rows, slong cols)
{
    slong i;

    if (rows != 0)
        mat->rows = flint_malloc(rows * sizeof(fmpq *));
    else
        mat->rows = NULL;

    mat->r = rows;
    mat->c = cols;

    if (rows != 0 && cols != 0)
    {
        slong num;
        int of;

        of = z_mul_checked(&num, rows, cols);

        if (of)
            flint_throw(FLINT_ERROR, "Overflow creating a %wd x %wd object\n", rows, cols);

        mat->entries = flint_calloc(num, sizeof(fmpq));

        /* Set denominators */
        for (i = 0; i < rows * cols; i++)
            mat->entries[i].den = WORD(1);

        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
    {
        mat->entries = NULL;
        if (rows != 0)
        {
            for (i = 0; i < rows; i++)
                mat->rows[i] = NULL;
        }
    }
}
