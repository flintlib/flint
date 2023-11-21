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
arb_mat_is_exact(const arb_mat_t A)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            if (!mag_is_zero(arb_radref(arb_mat_entry(A, i, j))))
                return 0;

    return 1;
}
