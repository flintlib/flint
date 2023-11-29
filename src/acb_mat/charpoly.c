/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"

void _acb_mat_charpoly(acb_ptr cp, const acb_mat_t mat, slong prec)
{
    const slong n = mat->r;

    if (n == 0)
    {
        acb_one(cp);
    }
    else if (n == 1)
    {
        acb_neg(cp + 0, acb_mat_entry(mat, 0, 0));
        acb_one(cp + 1);
    }
    else
    {
        slong i, k, t;
        acb_ptr a, A, s;

        a = _acb_vec_init(n * n);
        A = a + (n - 1) * n;

        _acb_vec_zero(cp, n + 1);
        acb_neg(cp + 0, acb_mat_entry(mat, 0, 0));

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                acb_set(a + 0 * n + i, acb_mat_entry(mat, i, t));
            }

            acb_set(A + 0, acb_mat_entry(mat, t, t));

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    acb_dot(s, NULL, 0, mat->rows[i], 1, a + (k - 1) * n, 1, t + 1, prec);
                }

                acb_set(A + k, a + k * n + t);
            }

            acb_dot(A + t, NULL, 0, mat->rows[t], 1, a + (t - 1) * n, 1, t + 1, prec);

            for (k = 0; k <= t; k++)
            {
                acb_dot(cp + k, cp + k, 1, A, 1, cp + k - 1, -1, k, prec);
                acb_sub(cp + k, cp + k, A + k, prec);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            acb_swap(cp + i, cp + (i - 1));

        acb_one(cp + 0);
        _acb_poly_reverse(cp, cp, n + 1, n + 1);
        _acb_vec_clear(a, n * n);
    }
}

void acb_mat_charpoly(acb_poly_t cp, const acb_mat_t mat, slong prec)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (acb_mat_charpoly).  Non-square matrix.\n");
    }

    acb_poly_fit_length(cp, mat->r + 1);
    _acb_poly_set_length(cp, mat->r + 1);
    _acb_mat_charpoly(cp->coeffs, mat, prec);
}
