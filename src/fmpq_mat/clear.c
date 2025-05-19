/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_clear(fmpq_mat_t mat)
{
    if (mat->entries != NULL)
    {
        slong i, j;

        for (i = 0; i < mat->r; i++)
            for (j = 0; j < mat->c; j++)
                fmpq_clear(fmpq_mat_entry(mat, i, j));

        flint_free(mat->entries);
    }
}
