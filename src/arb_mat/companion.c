/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_mat.h"

void
_arb_mat_companion(arb_mat_t A, arb_srcptr poly, slong prec)
{
    slong i, j, n;
    arb_t c;

    n = arb_mat_nrows(A);

    if (n == 0)
        return;

    for (i = 0; i < n - 1; i++)
        for (j = 0; j < n; j++)
            arb_set_ui(arb_mat_entry(A, i, j), (i + 1) == j);

    arb_init(c);
    arb_inv(c, poly + n, prec);
    arb_neg(c, c);
    for (j = 0; j < n; j++)
        arb_mul(arb_mat_entry(A, n - 1, j), poly + j, c, prec);
    arb_clear(c);
}

void
arb_mat_companion(arb_mat_t A, const arb_poly_t poly, slong prec)
{
    slong n = arb_mat_nrows(A);

    if (n != arb_poly_degree(poly) || n != arb_mat_ncols(A))
    {
        flint_throw(FLINT_ERROR, "arb_mat_companion: incompatible dimensions!\n");
    }

    _arb_mat_companion(A, poly->coeffs, prec);
}
