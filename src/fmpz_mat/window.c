/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_window_init(fmpz_mat_t window, const fmpz_mat_t mat, slong r1,
                     slong c1, slong r2, slong c2)
{
    slong i;
    window->entries = NULL;

    if (r2 > r1)
        window->rows = (fmpz **) flint_malloc((r2 - r1) * sizeof(fmpz *));
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
fmpz_mat_window_clear(fmpz_mat_t window)
{
    if (window->r != 0)
        flint_free(window->rows);
}

void
_fmpz_mat_window_readonly_init_strip_initial_zero_rows(fmpz_mat_t A, const fmpz_mat_t B)
{
    slong r = B->r;
    slong c = B->c;
    slong i;

    for (i = 0; i < r; i++)
    {
        if (!_fmpz_vec_is_zero(B->rows[i], c))
            break;
    }

    A->entries = NULL;
    A->rows = B->rows + i;
    A->r = r - i;
    A->c = c;
}
