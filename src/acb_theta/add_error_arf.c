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
arb_mat_add_error_arf(arb_mat_t mat, const arf_t err)
{
    slong k = acb_mat_nrows(mat);
    slong n = acb_mat_ncols(mat);
    slong i, j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < n; j++)
        {
            arb_add_error_arf(arb_mat_entry(mat, i, j), err);
        }
    }
}
