/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_inf_norm(arb_t res, const arb_mat_t A, slong prec)
{
    arb_t abs, sum;
    slong i, j;

    arb_init(abs);
    arb_init(sum);

    arb_zero(res);
    for (i = 0; i < arb_mat_nrows(A); i++)
    {
        arb_zero(sum);
        for (j = 0; j < arb_mat_ncols(A); j++)
        {
            arb_abs(abs, arb_mat_entry(A, i, j));
            arb_add(sum, sum, abs, prec);
        }
        arb_max(res, res, sum, prec);
    }

    arb_clear(abs);
    arb_clear(sum);
}
