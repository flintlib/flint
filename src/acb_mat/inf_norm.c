/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_inf_norm(arb_t res, const acb_mat_t A, slong prec)
{
    slong i, j, m, n;
    arb_t s, t;

    m = acb_mat_nrows(A);
    n = acb_mat_nrows(A);

    if (m == 0 || n == 0)
    {
        arb_zero(res);
        return;
    }

    arb_init(s);
    arb_init(t);

    arb_zero(res);

    for (i = 0; i < m; i++)
    {
        acb_abs(s, acb_mat_entry(A, i, 0), prec);

        for (j = 1; j < n; j++)
        {
            acb_abs(t, acb_mat_entry(A, i, j), prec);
            arb_add(s, s, t, prec);
        }

        arb_max(res, res, s, prec);
    }

    arb_clear(s);
    arb_clear(t);
}
