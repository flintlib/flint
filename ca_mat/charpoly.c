/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
_ca_mat_charpoly(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    const slong n = mat->r;

    if (n == 0)
    {
        ca_one(cp, ctx);
    }
    else if (n == 1)
    {
        ca_neg(cp + 0, ca_mat_entry(mat, 0, 0), ctx);
        ca_one(cp + 1, ctx);
    }
    else
    {
        slong i, k, t;
        ca_ptr a, A, s;

        a = _ca_vec_init(n * n, ctx);
        A = a + (n - 1) * n;

        _ca_vec_zero(cp, n + 1, ctx);
        ca_neg(cp + 0, ca_mat_entry(mat, 0, 0), ctx);

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                ca_set(a + 0 * n + i, ca_mat_entry(mat, i, t), ctx);
            }

            ca_set(A + 0, ca_mat_entry(mat, t, t), ctx);

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    ca_dot(s, NULL, 0, mat->rows[i], 1, a + (k - 1) * n, 1, t + 1, ctx);
                }

                ca_set(A + k, a + k * n + t, ctx);
            }

            ca_dot(A + t, NULL, 0, mat->rows[t], 1, a + (t - 1) * n, 1, t + 1, ctx);

            for (k = 0; k <= t; k++)
            {
                ca_dot(cp + k, cp + k, 1, A, 1, cp + k - 1, -1, k, ctx);
                ca_sub(cp + k, cp + k, A + k, ctx);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            ca_swap(cp + i, cp + (i - 1), ctx);

        ca_one(cp + 0, ctx);
        _ca_poly_reverse(cp, cp, n + 1, n + 1, ctx);
        _ca_vec_clear(a, n * n, ctx);
    }
}

void ca_mat_charpoly(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    if (mat->r != mat->c)
    {
        flint_printf("Exception (ca_mat_charpoly).  Non-square matrix.\n");
        flint_abort();
    }

    ca_poly_fit_length(cp, mat->r + 1, ctx);
    _ca_poly_set_length(cp, mat->r + 1, ctx);
    _ca_mat_charpoly(cp->coeffs, mat, ctx);
}
