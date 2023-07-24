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
arb_mat_bilinear_form(arb_t x, const arb_mat_t A, arb_srcptr v1, arb_srcptr v2, slong prec)
{
    slong nrow = arb_mat_nrows(A);
    slong ncol = arb_mat_ncols(A);
    arb_mat_t col, row, prod, scal;
    slong k;

    arb_mat_init(col, ncol, 1);
    arb_mat_init(row, 1, nrow);
    arb_mat_init(prod, nrow, 1);
    arb_mat_init(scal, 1, 1);

    for (k = 0; k < nrow; k++)
    {
        arb_set(arb_mat_entry(row, 0, k), &v1[k]);
    }
    for (k = 0; k < ncol; k++)
    {
        arb_set(arb_mat_entry(col, k, 0), &v2[k]);
    }
    arb_mat_mul(prod, A, col, prec);
    arb_mat_mul(scal, row, prod, prec);
    arb_set(x, arb_mat_entry(scal, 0, 0));

    arb_mat_clear(col);
    arb_mat_clear(row);
    arb_mat_clear(prod);
    arb_mat_clear(scal);
}
