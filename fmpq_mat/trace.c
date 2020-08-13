/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void
fmpq_mat_trace(fmpq_t trace, const fmpq_mat_t mat)
{
    slong i, n = fmpq_mat_nrows(mat);

    if (n == 0)
        fmpq_zero(trace);
    else
    {
        fmpq_set(trace, fmpq_mat_entry(mat, 0, 0));
        for (i = 1; i < n; i++)
            fmpq_add(trace, trace, fmpq_mat_entry(mat, i, i));
    }
}
