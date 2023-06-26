/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_mat_get_real(arb_mat_t re, const acb_mat_t mat)
{
    slong nrows = acb_mat_nrows(mat);
    slong ncols = acb_mat_ncols(mat);
    slong i, j;

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            acb_get_real(arb_mat_entry(re, i, j), acb_mat_entry(mat, i, j));
        }
    }
}
