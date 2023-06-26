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
acb_mat_ninf(arb_t norm, const acb_mat_t mat, slong prec)
{
    slong n = acb_mat_nrows(mat);
    slong k = acb_mat_ncols(mat);
    arb_t abs, sum;
    slong i, j;

    arb_init(abs);
    arb_init(sum);

    arb_zero(norm);
    for (i = 0; i < n; i++)
    {
        arb_zero(sum);
        for (j = 0; j < k; j++)
        {
            acb_abs(abs, acb_mat_entry(mat, i, j), prec);
            arb_add(sum, sum, abs, prec);
        }
        arb_max(norm, norm, sum, prec);
    }

    arb_clear(abs);
    arb_clear(sum);
}
