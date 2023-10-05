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
acb_mat_vector_mul_row(acb_ptr res, acb_srcptr v, const acb_mat_t A, slong prec)
{
    slong nrow = acb_mat_nrows(A);
    slong ncol = acb_mat_ncols(A);
    acb_mat_t r, p;
    slong k;

    acb_mat_init(r, 1, nrow);
    acb_mat_init(p, 1, ncol);

    for (k = 0; k < nrow; k++)
    {
        acb_set(acb_mat_entry(r, 0, k), &v[k]);
    }
    acb_mat_mul(p, r, A, prec);
    for (k = 0; k < ncol; k++)
    {
        acb_set(&res[k], acb_mat_entry(p, 0, k));
    }

    acb_mat_clear(r);
    acb_mat_clear(p);
}

void
acb_mat_vector_mul_col(acb_ptr res, const acb_mat_t A, acb_srcptr v, slong prec)
{
    slong nrow = acb_mat_nrows(A);
    slong ncol = acb_mat_ncols(A);
    acb_mat_t c, p;
    slong k;

    acb_mat_init(c, ncol, 1);
    acb_mat_init(p, nrow, 1);

    for (k = 0; k < ncol; k++)
    {
        acb_set(acb_mat_entry(c, k, 0), &v[k]);
    }
    acb_mat_mul(p, A, c, prec);
    for (k = 0; k < nrow; k++)
    {
        acb_set(&res[k], acb_mat_entry(p, k, 0));
    }

    acb_mat_clear(c);
    acb_mat_clear(p);
}
