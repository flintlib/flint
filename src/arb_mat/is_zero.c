/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
arb_mat_is_zero(const arb_mat_t A)
{
    slong i, j, n, m;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(A);

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            if (!arb_is_zero(arb_mat_entry(A, i, j)))
                return 0;

    return 1;
}
