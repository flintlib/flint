/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
_acb_mat_diag_prod(acb_t res, const acb_mat_t A, slong a, slong b, slong prec)
{
    if (b - a == 0)
        acb_one(res);
    else if (b - a == 1)
        acb_set_round(res, acb_mat_entry(A, a, a), prec);
    else if (b - a == 2)
        acb_mul(res, acb_mat_entry(A, a, a), acb_mat_entry(A, a + 1, a + 1), prec);
    else if (b - a == 3)
    {
        acb_mul(res, acb_mat_entry(A, a, a), acb_mat_entry(A, a + 1, a + 1), prec);
        acb_mul(res, res, acb_mat_entry(A, a + 2, a + 2), prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        _acb_mat_diag_prod(t, A, a, a + (b - a) / 2, prec);
        _acb_mat_diag_prod(res, A, a + (b - a) / 2, b, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }
}

void
acb_mat_diag_prod(acb_t res, const acb_mat_t A, slong prec)
{
    slong m, n;

    m = acb_mat_nrows(A);
    n = acb_mat_nrows(A);

    _acb_mat_diag_prod(res, A, 0, FLINT_MIN(m, n), prec);
}
