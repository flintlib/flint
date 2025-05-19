/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
_fmpz_mat_window_readonly_init_strip_initial_zero_rows(fmpz_mat_t A, const fmpz_mat_t B)
{
    slong r = B->r;
    slong c = B->c;
    slong i;

    for (i = 0; i < r; i++)
    {
        if (!_fmpz_vec_is_zero(fmpz_mat_row(B, i), c))
            break;
    }

    fmpz_mat_window_init(A, B, i, 0, r, c);
}
