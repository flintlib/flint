/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_dirichlet.h"

void
arb_mat_dct(arb_mat_t res, int kind, slong prec)
{
    acb_dirichlet_roots_t roots;
    acb_t t;
    arb_t v;
    slong n, r, c, i, j;

    r = arb_mat_nrows(res);
    c = arb_mat_ncols(res);
    n = FLINT_MIN(r, c);

    if (n == 0)
        return;

    acb_dirichlet_roots_init(roots, 4 * n, (r - 1) * c, prec);
    acb_init(t);
    arb_init(v);

    arb_set_ui(v, n);
    arb_rsqrt(v, v, prec);

    if (r > 0)
    {
        for (j = 0; j < c; j++)
            arb_set(arb_mat_entry(res, 0, j), v);
    }

    arb_set_ui(v, n);
    arb_mul_2exp_si(v, v, -1);
    arb_rsqrt(v, v, prec);

    for (i = 1; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            acb_dirichlet_root(t, roots, i * (2 * j + 1), prec);
            arb_mul(arb_mat_entry(res, i, j), acb_realref(t), v, prec);
        }
    }

    acb_dirichlet_roots_clear(roots);
    acb_clear(t);
    arb_clear(v);
}

