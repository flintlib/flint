/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_find_pivot_any(const fmpz_mat_t mat,
                                    slong start_row, slong end_row, slong c)
{
    slong r;

    for (r = start_row; r < end_row; r++)
    {
        if (!fmpz_is_zero(fmpz_mat_entry(mat, r, c)))
            return r;
    }

    return -1;
}
