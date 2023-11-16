/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_frobenius_norm(arb_t res, const arb_mat_t A, slong prec)
{
    slong i, j, r, c;

    r = arb_mat_nrows(A);
    c = arb_mat_ncols(A);

    arb_zero(res);

    if (r == 0 || c == 0)
        return;

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            arb_srcptr x = arb_mat_entry(A, i, j);
            arb_addmul(res, x, x, prec);
        }
    }

    arb_sqrtpos(res, res, prec);
}
