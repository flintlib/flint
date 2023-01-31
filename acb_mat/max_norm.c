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
acb_mat_max_norm(arb_t res, const acb_mat_t A, slong prec)
{
    slong i, j;
    arb_t abs;

    arb_init(abs);
    arb_zero(res);
    for (i = 0; i < acb_mat_nrows(A); i++)
    {
        for (j = 0; j < acb_mat_ncols(A); j++)
        {
            acb_abs(abs, acb_mat_entry(mat, i, j));
            arb_max(res, res, abs, prec);
        }
    }
    arb_clear(abs);
}
