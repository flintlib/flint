/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_clear(fmpq_mat_t mat)
{
    if (mat->entries)
    {
        slong i;

        for (i = 0; i < mat->r * mat->c; i++)
            fmpq_clear(mat->entries + i);

        flint_free(mat->entries);
        flint_free(mat->rows);
    } else if (mat->r != 0)
        flint_free(mat->rows);
}
