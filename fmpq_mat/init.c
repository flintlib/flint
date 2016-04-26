/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_init(fmpq_mat_t mat, slong rows, slong cols)
{
    if ((rows) && (cols))
    {
        slong i;
        mat->entries = (fmpq *) flint_calloc(rows * cols, sizeof(fmpq));
        mat->rows = (fmpq **) flint_malloc(rows * sizeof(fmpq *));

        /* Set denominators */
        for (i = 0; i < rows * cols; i++)
            mat->entries[i].den = WORD(1);

        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
        mat->entries = NULL;

    mat->r = rows;
    mat->c = cols;
}
