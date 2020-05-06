/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpq_mat.h"

void
fmpq_mat_window_init(fmpq_mat_t window, const fmpq_mat_t mat, slong r1,
                     slong c1, slong r2, slong c2)
{
    slong i;
    window->entries = NULL;

    if (r2 > r1 && c2 > c1)
    {
        window->rows = (fmpq **) flint_malloc((r2 - r1) * sizeof(fmpq *));
        window->r = r2 - r1;
        window->c = c2 - c1;

        if (mat->c > 0)
        {
            for (i = 0; i < r2 - r1; i++)
                window->rows[i] = mat->rows[r1 + i] + c1;
        }
    }
    else
        memset(window, 0, sizeof(*window));
}
