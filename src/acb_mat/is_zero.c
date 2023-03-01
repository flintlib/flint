/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_is_zero(const acb_mat_t A)
{
    slong i, j, n, m;

    n = acb_mat_nrows(A);
    m = acb_mat_ncols(A);

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            if (!acb_is_zero(acb_mat_entry(A, i, j)))
                return 0;

    return 1;
}
