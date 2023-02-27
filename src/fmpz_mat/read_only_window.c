/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void _fmpz_mat_read_only_window_init_strip_initial_zero_rows(
    fmpz_mat_t A,
    const fmpz_mat_t B)
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
