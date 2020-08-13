/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

/*
    Assumes that \code{mat} is an $n \times n$ matrix and sets \code{(cp,n+1)} 
    to its characteristic polynomial.

    Employs a division-free algorithm using $O(n^4)$ ring operations.
 */

void _fmpz_mat_charpoly_berkowitz(fmpz *cp, const fmpz_mat_t mat)
{
    const slong n = mat->r;

    if (n == 0)
    {
        fmpz_one(cp);
    }
    else if (n == 1)
    {
        fmpz_neg(cp + 0, fmpz_mat_entry(mat, 0, 0));
        fmpz_one(cp + 1);
    }
    else
    {
        slong i, j, k, t;
        fmpz *a, *A, *s;

        a = _fmpz_vec_init(n * n);
        A = a + (n - 1) * n;

        _fmpz_vec_zero(cp, n + 1);
        fmpz_neg(cp + 0, fmpz_mat_entry(mat, 0, 0));

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                fmpz_set(a + 0 * n + i, fmpz_mat_entry(mat, i, t));
            }

            fmpz_set(A + 0, fmpz_mat_entry(mat, t, t));

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = a + k * n + i;
                    fmpz_zero(s);
                    for (j = 0; j <= t; j++)
                    {
                        fmpz_addmul(s, fmpz_mat_entry(mat, i, j), a + (k - 1) * n + j);
                    }
                }
                fmpz_set(A + k, a + k * n + t);
            }

            fmpz_zero(A + t);
            for (j = 0; j <= t; j++)
            {
                fmpz_addmul(A + t, fmpz_mat_entry(mat, t, j), a + (t - 1) * n + j);
            }

            for (k = 0; k <= t; k++)
            {
                for (j = 0; j < k; j++)
                {
                    fmpz_submul(cp + k, A + j, cp + (k - j - 1));
                }
                fmpz_sub(cp + k, cp + k, A + k);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
        {
            fmpz_swap(cp + i, cp + (i - 1));
        }
        fmpz_one(cp + 0);

        _fmpz_poly_reverse(cp, cp, n + 1, n + 1);

        _fmpz_vec_clear(a, n * n);
    }
}

void fmpz_mat_charpoly_berkowitz(fmpz_poly_t cp, const fmpz_mat_t mat)
{
     fmpz_poly_fit_length(cp, mat->r + 1);
    _fmpz_poly_set_length(cp, mat->r + 1);

    _fmpz_mat_charpoly_berkowitz(cp->coeffs, mat);
}

