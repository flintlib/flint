/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_init(fmpz_mat_t mat, slong rows, slong cols)
{
    slong i;
    
    if (rows != 0)
        mat->rows = (fmpz **) flint_malloc(rows * sizeof(fmpz *));
    else
        mat->rows = NULL;
    
    if (rows != 0 && cols != 0)       /* Allocate space for r*c small entries */
    {
        mat->entries = (fmpz *) flint_calloc(flint_mul_sizes(rows, cols), sizeof(fmpz));

        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
    {
        mat->entries = NULL;
        for (i = 0; i < rows; i++)
            mat->rows[i] = NULL;
    }

    mat->r = rows;
    mat->c = cols;
}
