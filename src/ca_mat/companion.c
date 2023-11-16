/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
_ca_mat_companion(ca_mat_t A, ca_srcptr poly, const ca_t c, ca_ctx_t ctx)
{
    slong i, j, n;

    n = ca_mat_nrows(A);

    if (n == 0)
        return;

    for (i = 0; i < n - 1; i++)
        for (j = 0; j < n; j++)
            ca_set_ui(ca_mat_entry(A, i, j), (i + 1) == j, ctx);

    for (j = 0; j < n; j++)
        ca_mul(ca_mat_entry(A, n - 1, j), poly + j, c, ctx);
}

int
ca_mat_companion(ca_mat_t A, const ca_poly_t poly, ca_ctx_t ctx)
{
    ca_t c;
    int res;
    slong n = ca_mat_nrows(A);

    if (n != poly->length - 1 || n != ca_mat_ncols(A))
    {
        return 0;
    }

    if (CA_IS_SPECIAL(poly->coeffs + n))
        return 0;

    ca_init(c, ctx);

    ca_inv(c, poly->coeffs + n, ctx);
    ca_neg(c, c, ctx);

    if (CA_IS_SPECIAL(c))
    {
        res = 0;
    }
    else
    {
        _ca_mat_companion(A, poly->coeffs, c, ctx);
        res = 1;
    }

    ca_clear(c, ctx);

    return res;
}
