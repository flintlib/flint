/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_dirichlet.h"

void
acb_mat_dft(acb_mat_t res, int kind, slong prec)
{
    acb_dirichlet_roots_t roots;
    acb_t t;
    acb_t v;
    slong n, r, c, i, j;

    r = arb_mat_nrows(res);
    c = arb_mat_ncols(res);
    n = FLINT_MIN(r, c);

    if (n == 0)
        return;

    acb_dirichlet_roots_init(roots, n, (r - 1) * c, prec);
    acb_init(t);
    acb_init(v);

    acb_set_ui(v, n);
    acb_rsqrt(v, v, prec);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            acb_dirichlet_root(t, roots, i * j, prec);
            acb_conj(t, t);
            acb_mul(acb_mat_entry(res, i, j), t, v, prec);
        }
    }

    acb_dirichlet_roots_clear(roots);
    acb_clear(t);
    acb_clear(v);
}

