/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2015 Elena Sergeichave

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_window_init(fmpz_poly_mat_t window, const fmpz_poly_mat_t mat, slong r1,
                     slong c1, slong r2, slong c2)
{
    slong i;
    window->entries = NULL;

    if (r2 > r1)
        window->rows = (fmpz_poly_struct **) flint_malloc((r2 - r1)
                                                  * sizeof(fmpz_poly_struct *));
    else
        window->rows = NULL;

    if (mat->c > 0)
    {
        for (i = 0; i < r2 - r1; i++)
            window->rows[i] = mat->rows[r1 + i] + c1;
    } else
    {
        for (i = 0; i < r2 - r1; i++)
            window->rows[i] = NULL;
    }

    window->r = r2 - r1;
    window->c = c2 - c1;
}

void
fmpz_poly_mat_window_clear(fmpz_poly_mat_t window)
{
    if (window->r != 0)
        flint_free(window->rows);
}
