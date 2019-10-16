/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

void _fmpz_mod_mat_reduce(fmpz_mod_mat_t mat)
{
    slong i, j;
    fmpz *entry;
    for (i = 0; i < mat->mat->r; i++)
    {
        for (j = 0; j < mat->mat->c; j++)
        {
            entry = fmpz_mod_mat_entry(mat, i, j);
            fmpz_mod(entry, entry, mat->mod);
        }
    }
}
