/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "ca_mat.h"

void
_ca_mat_ca_poly_evaluate(ca_mat_t y, ca_srcptr poly,
    slong len, const ca_mat_t x, ca_ctx_t ctx)
{
    slong i, j, m, r, n;
    ca_mat_struct * xs;
    ca_mat_t s, t;

    if (len == 0)
    {
        ca_mat_zero(y, ctx);
        return;
    }

    if (len == 1)
    {
        ca_mat_set_ca(y, poly + 0, ctx);
        return;
    }

    if (len == 2)
    {
        ca_mat_mul_ca(y, x, poly + 1, ctx);
        ca_mat_add_ca(y, y, poly + 0, ctx);
        return;
    }

    n = ca_mat_nrows(x);

    m = n_sqrt(len) + 1;
    r = (len + m - 1) / m;

    xs = flint_malloc(sizeof(ca_mat_struct) * (m + 1));
    for (i = 0; i <= m; i++)
    {
        ca_mat_init(xs + i, n, n, ctx);

        if (i == 0)
            ca_mat_one(xs + i, ctx);
        else if (i == 1)
            ca_mat_set(xs + i, x, ctx);
        else
            ca_mat_mul(xs + i, xs + i - 1, x, ctx);
    }

    ca_mat_init(s, n, n, ctx);
    ca_mat_init(t, n, n, ctx);

    ca_mat_set_ca(y, poly + (r - 1) * m, ctx);
    for (j = 1; (r - 1) * m + j < len; j++)
        ca_mat_addmul_ca(y, xs + j, poly + (r - 1) * m + j, ctx);

    for (i = r - 2; i >= 0; i--)
    {
        ca_mat_set_ca(s, poly + i * m, ctx);
        for (j = 1; j < m; j++)
            ca_mat_addmul_ca(s, xs + j, poly + i * m + j, ctx);

        ca_mat_mul(y, y, xs + m, ctx);
        ca_mat_add(y, y, s, ctx);
    }

    for (i = 0; i <= m; i++)
        ca_mat_clear(xs + i, ctx);
    flint_free(xs);
    ca_mat_clear(s, ctx);
    ca_mat_clear(t, ctx);
}

void
ca_mat_ca_poly_evaluate(ca_mat_t res, const ca_poly_t f, const ca_mat_t a, ca_ctx_t ctx)
{
    _ca_mat_ca_poly_evaluate(res, f->coeffs, f->length, a, ctx);
}
