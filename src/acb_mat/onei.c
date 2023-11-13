/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_onei(acb_mat_t mat)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(mat); i++)
    {
        for (j = 0; j < acb_mat_ncols(mat); j++)
        {
            if (i == j)
            {
                acb_onei(acb_mat_entry(mat, i, j));
            }
            else
            {
                acb_zero(acb_mat_entry(mat, i, j));
            }
        }
    }
}
