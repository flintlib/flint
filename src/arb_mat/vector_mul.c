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
arb_mat_vector_mul_row(arb_ptr res, arb_srcptr v, const arb_mat_t A, slong prec)
{
    slong nrow = arb_mat_nrows(A);
    slong ncol = arb_mat_ncols(A);
    arb_mat_t r, p;
    slong k;

    arb_mat_init(r, 1, nrow);
    arb_mat_init(p, 1, ncol);

    for (k = 0; k < nrow; k++)
    {
        arb_set(arb_mat_entry(r, 0, k), &v[k]);
    }
    arb_mat_mul(p, r, A, prec);
    for (k = 0; k < ncol; k++)
    {
        arb_set(&res[k], arb_mat_entry(p, 0, k));
    }

    arb_mat_clear(r);
    arb_mat_clear(p);
}

void
arb_mat_vector_mul_col(arb_ptr res, const arb_mat_t A, arb_srcptr v, slong prec)
{
    slong nrow = arb_mat_nrows(A);
    slong ncol = arb_mat_ncols(A);
    arb_mat_t c, p;
    slong k;

    arb_mat_init(c, ncol, 1);
    arb_mat_init(p, nrow, 1);

    for (k = 0; k < ncol; k++)
    {
        arb_set(arb_mat_entry(c, k, 0), &v[k]);
    }
    arb_mat_mul(p, A, c, prec);
    for (k = 0; k < nrow; k++)
    {
        arb_set(&res[k], arb_mat_entry(p, k, 0));
    }

    arb_mat_clear(c);
    arb_mat_clear(p);
}
