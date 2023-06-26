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
acb_mat_set_arb_arb(acb_mat_t mat, const arb_mat_t re, const arb_mat_t im)
{
    slong nrows = acb_mat_nrows(re);
    slong ncols = acb_mat_ncols(re);
    slong i, j;

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            acb_set_arb_arb(acb_mat_entry(mat, i, j),
                            arb_mat_entry(re, i, j), arb_mat_entry(im, i, j));
        }
    }
}
