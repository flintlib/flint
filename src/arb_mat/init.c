/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_init(arb_mat_t mat, slong r, slong c)
{
    if (r != 0 && c != 0)
    {
        slong i;

        mat->entries = _arb_vec_init(r * c);
        mat->rows = (arb_ptr *) flint_malloc(r * sizeof(arb_ptr));

        for (i = 0; i < r; i++)
            mat->rows[i] = mat->entries + i * c;
    }
    else
        mat->entries = NULL;

    mat->r = r;
    mat->c = c;
}
