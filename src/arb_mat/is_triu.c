/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
arb_mat_is_triu(const arb_mat_t A)
{
    slong i, j, n, m;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(A);

    for (i = 1; i < n; i++)
        for (j = 0; j < FLINT_MIN(i, m); j++)
            if (!arb_is_zero(arb_mat_entry(A, i, j)))
                return 0;

    return 1;
}
