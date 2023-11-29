/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"

void
_acb_mat_companion(acb_mat_t A, acb_srcptr poly, slong prec)
{
    slong i, j, n;
    acb_t c;

    n = acb_mat_nrows(A);

    if (n == 0)
        return;

    for (i = 0; i < n - 1; i++)
        for (j = 0; j < n; j++)
            acb_set_ui(acb_mat_entry(A, i, j), (i + 1) == j);

    acb_init(c);
    acb_inv(c, poly + n, prec);
    acb_neg(c, c);
    for (j = 0; j < n; j++)
        acb_mul(acb_mat_entry(A, n - 1, j), poly + j, c, prec);
    acb_clear(c);
}

void
acb_mat_companion(acb_mat_t A, const acb_poly_t poly, slong prec)
{
    slong n = acb_mat_nrows(A);

    if (n != acb_poly_degree(poly) || n != acb_mat_ncols(A))
    {
        flint_throw(FLINT_ERROR, "acb_mat_companion: incompatible dimensions!\n");
    }

    _acb_mat_companion(A, poly->coeffs, prec);
}
