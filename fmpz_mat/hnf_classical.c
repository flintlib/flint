/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_hnf_classical(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong i, i0, j, j2, k, l, m, n;
    fmpz_t min, q;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_init(q);
    fmpz_mat_set(H, A);

    for (j = 0, k = 0, l = (n - m) * (n > m); n - j != l; j++, k++)
    {
        int col_finished = 1;
        for (i = k + 1; (i < m) && col_finished; i++)
            col_finished = fmpz_is_zero(fmpz_mat_entry(H, i, j));
        if (col_finished)
        {
            if (fmpz_sgn(fmpz_mat_entry(H, k, j)) < 0)
            {
                for (j2 = j; j2 < n; j2++)
                    fmpz_neg(fmpz_mat_entry(H, k, j2),
                             fmpz_mat_entry(H, k, j2));
            }
            if (fmpz_is_zero(fmpz_mat_entry(H, k, j)))
            {
                k--;
                if (l > 0)
                    l--;
            }
            else
            {
                /* reduce first entries of column j with row k */
                for (i = 0; i < k; i++)
                {
                    fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                                fmpz_mat_entry(H, k, j));
                    for (j2 = j; j2 < n; j2++)
                    {
                        fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                                    fmpz_mat_entry(H, k, j2));
                    }
                }
            }
        }
        else
        {
            i0 = 0;
            fmpz_init(min);
            /* locate non-zero entry in column j below k with lowest absolute
               value */
            for (i = k + 1; i < m; i++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(H, i, j)))
                    continue;
                if (fmpz_is_zero(min) ||
                    fmpz_cmpabs(fmpz_mat_entry(H, i, j), min) < 0)
                {
                    i0 = i;
                    fmpz_abs(min, fmpz_mat_entry(H, i, j));
                }
            }
            /* move the row found to row k */
            if (i0 > k)
                fmpz_mat_swap_rows(H, NULL, i0, k);

            if (fmpz_sgn(fmpz_mat_entry(H, k, j)) < 0)
            {
                for (j2 = j; j2 < n; j2++)
                {
                    fmpz_neg(fmpz_mat_entry(H, k, j2),
                             fmpz_mat_entry(H, k, j2));
                }
            }
            /* reduce lower entries of column j with row k */
            for (i = k + 1; i < m; i++)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                            fmpz_mat_entry(H, k, j));
                for (j2 = j; j2 < n; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                                fmpz_mat_entry(H, k, j2));
                }
            }
            /* don't move to the next column yet */
            j--;
            k--;
            fmpz_clear(min);
        }
    }

    fmpz_clear(q);
}
