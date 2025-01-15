/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_clear(fmpz_mat_t mat)
{
    if (mat->entries != NULL)
    {
        slong i;
        for (i = 0; i < mat->r * mat->c; i++)
            fmpz_clear(mat->entries + i);   /* Clear all coefficients */
        flint_free(mat->entries);
    }
}
