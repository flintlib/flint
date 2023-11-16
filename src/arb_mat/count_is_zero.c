/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

slong
arb_mat_count_is_zero(const arb_mat_t mat)
{
    slong nz;
    slong i, j;
    nz = 0;
    for (i = 0; i < arb_mat_nrows(mat); i++)
    {
        for (j = 0; j < arb_mat_ncols(mat); j++)
        {
            if (arb_is_zero(arb_mat_entry(mat, i, j)))
            {
                nz++;
            }
        }
    }
    return nz;
}
