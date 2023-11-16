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
ca_mat_adjugate_charpoly(ca_mat_t adj, ca_t det, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_poly_t pol;
    slong n;

    n = ca_mat_nrows(A);

    if (n == 0)
    {
        ca_one(det, ctx);
        return;
    }

    ca_poly_init(pol, ctx);
    ca_mat_charpoly(pol, A, ctx);
    ca_swap(det, ca_poly_coeff_ptr(pol, 0), ctx);
    ca_poly_shift_right(pol, pol, 1, ctx);
    ca_mat_ca_poly_evaluate(adj, pol, A, ctx);

    if (n % 2)
        ca_neg(det, det, ctx);
    else
        ca_mat_neg(adj, adj, ctx);

    ca_poly_clear(pol, ctx);
}
