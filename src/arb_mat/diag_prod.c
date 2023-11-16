/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
_arb_mat_diag_prod(arb_t res, const arb_mat_t A, slong a, slong b, slong prec)
{
    if (b - a == 0)
    {
        arb_one(res);
    }
    else if (b - a == 1)
    {
        arb_set_round(res, arb_mat_entry(A, a, a), prec);
    }
    else
    {
        slong i;

        arb_mul(res, arb_mat_entry(A, a, a), arb_mat_entry(A, a + 1, a + 1), prec);
        for (i = a + 2; i < b; i++)
            arb_mul(res, res, arb_mat_entry(A, i, i), prec);

        /* no advantage? */

        /*
        arb_t t;
        arb_init(t);
        _arb_mat_diag_prod(t, A, a, a + (b - a) / 2, prec);
        _arb_mat_diag_prod(res, A, a + (b - a) / 2, b, prec);
        arb_mul(res, res, t, prec);
        arb_clear(t);
        */
    }
}

void
arb_mat_diag_prod(arb_t res, const arb_mat_t A, slong prec)
{
    slong m, n;

    m = arb_mat_nrows(A);
    n = arb_mat_nrows(A);

    _arb_mat_diag_prod(res, A, 0, FLINT_MIN(m, n), prec);
}
