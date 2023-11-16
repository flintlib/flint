/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_frobenius_norm(arb_t res, const acb_mat_t A, slong prec)
{
    slong i, j, r, c;

    r = acb_mat_nrows(A);
    c = acb_mat_ncols(A);

    arb_zero(res);

    if (r == 0 || c == 0)
        return;

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            acb_srcptr z = acb_mat_entry(A, i, j);
            arb_addmul(res, acb_realref(z), acb_realref(z), prec);
            arb_addmul(res, acb_imagref(z), acb_imagref(z), prec);
        }
    }

    arb_sqrtpos(res, res, prec);
}
