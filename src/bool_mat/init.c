/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

void
bool_mat_init(bool_mat_t mat, slong r, slong c)
{
    mat->entries = NULL;
    mat->rows = NULL;
    mat->r = r;
    mat->c = c;
    if (r != 0 && c != 0)
    {
        slong i;
        mat->entries = flint_calloc(r * c, sizeof(int));
        mat->rows = flint_malloc(r * sizeof(int *));
        for (i = 0; i < r; i++)
            mat->rows[i] = mat->entries + i * c;
    }
}
