/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_mat.h"

void _arb_mat_charpoly(arb_ptr cp, const arb_mat_t mat, slong prec)
{
    const slong n = mat->r;

    if (n == 0)
    {
        arb_one(cp);
    }
    else if (n == 1)
    {
        arb_neg(cp + 0, arb_mat_entry(mat, 0, 0));
        arb_one(cp + 1);
    }
    else
    {
        slong i, k, t;
        arb_ptr a, A, s;

        a = _arb_vec_init(n * n);
        A = a + (n - 1) * n;

        _arb_vec_zero(cp, n + 1);
        arb_neg(cp + 0, arb_mat_entry(mat, 0, 0));

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                arb_set(a + 0 * n + i, arb_mat_entry(mat, i, t));
            }

            arb_set(A + 0, arb_mat_entry(mat, t, t));

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    arb_dot(s, NULL, 0, mat->rows[i], 1, a + (k - 1) * n, 1, t + 1, prec);
                }

                arb_set(A + k, a + k * n + t);
            }

            arb_dot(A + t, NULL, 0, mat->rows[t], 1, a + (t - 1) * n, 1, t + 1, prec);

            for (k = 0; k <= t; k++)
            {
                arb_dot(cp + k, cp + k, 1, A, 1, cp + k - 1, -1, k, prec);
                arb_sub(cp + k, cp + k, A + k, prec);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            arb_swap(cp + i, cp + (i - 1));

        arb_one(cp + 0);
        _arb_poly_reverse(cp, cp, n + 1, n + 1);
        _arb_vec_clear(a, n * n);
    }
}

void arb_mat_charpoly(arb_poly_t cp, const arb_mat_t mat, slong prec)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (arb_mat_charpoly).  Non-square matrix.\n");
    }

    arb_poly_fit_length(cp, mat->r + 1);
    _arb_poly_set_length(cp, mat->r + 1);
    _arb_mat_charpoly(cp->coeffs, mat, prec);
}
