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
ca_mat_dft(ca_mat_t res, int type, ca_ctx_t ctx)
{
    ca_ptr w;
    slong n, r, c, i, j;

    r = ca_mat_nrows(res);
    c = ca_mat_ncols(res);
    n = FLINT_MIN(r, c);

    if (n == 0)
        return;

    w = _ca_vec_init(2 * n, ctx);

    for (i = 0; i < 2 * n; i++)
    {
        if (i == 0)
        {
            ca_one(w + i, ctx);
        }
        else if (i == 1)
        {
            ca_pi_i(w + i, ctx);
            ca_mul_ui(w + i, w + i, 2, ctx);
            ca_div_si(w + i, w + i, n, ctx);
            ca_exp(w + i, w + i, ctx);
            if (type == 0 || type == 2)
                ca_inv(w + i, w + i, ctx);
        }
        else
        {
            ca_mul(w + i, w + i - 1, w + 1, ctx);
        }
    }

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            ca_set(ca_mat_entry(res, i, j), w + (i * j) % (2 * n), ctx);
        }
    }

    _ca_vec_clear(w, 2 * n, ctx);

    if (type == 1)
    {
        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                ca_div_ui(ca_mat_entry(res, i, j), ca_mat_entry(res, i, j), n, ctx);
    }
    else if (type == 2 || type == 3)
    {
        ca_t t;
        ca_init(t, ctx);

        ca_sqrt_ui(t, n, ctx);
        ca_inv(t, t, ctx);

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                ca_mul(ca_mat_entry(res, i, j), ca_mat_entry(res, i, j), t, ctx);

        ca_clear(t, ctx);
    }
}
