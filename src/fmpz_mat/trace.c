/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_trace(fmpz_t trace, const fmpz_mat_t mat)
{
    slong i, n = fmpz_mat_nrows(mat);

    if (n == 0)
        fmpz_zero(trace);
    else
    {
        fmpz_set(trace, fmpz_mat_entry(mat, 0, 0));
        for (i = 1; i < n; i++)
            fmpz_add(trace, trace, fmpz_mat_entry(mat, i, i));
    }
}
