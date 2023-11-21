/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_set(acb_mat_t dest, const acb_mat_t src)
{
    slong i, j;

    if (dest != src && acb_mat_ncols(src) != 0)
    {
        for (i = 0; i < acb_mat_nrows(src); i++)
            for (j = 0; j < acb_mat_ncols(src); j++)
                acb_set(acb_mat_entry(dest, i, j),
                    acb_mat_entry(src, i, j));
    }
}
