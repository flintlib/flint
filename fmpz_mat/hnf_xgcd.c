/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void fmpz_mat_hnf_xgcd(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong j, j2, i, k, l;
    fmpz_t r1d, r2d, b, u, v, d, q;

    fmpz_init(r1d);
    fmpz_init(r2d);
    fmpz_init(b);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(d);
    fmpz_init(q);
    fmpz_mat_set(H, A);
    for (j = 0, k = 0, l = (A->c - A->r)*(A->c > A->r); A->c - j != l; j++, k++)
    {
        for (i = k + 1; i != A->r; i++)
        {
            /* reduce row i - 1 with row i */
            if (fmpz_is_zero(fmpz_mat_entry(H, i - 1, j)))
                continue;
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, i, j),
                    fmpz_mat_entry(H, i - 1, j));
            fmpz_divexact(r2d, fmpz_mat_entry(H, i - 1, j), d);
            fmpz_divexact(r1d, fmpz_mat_entry(H, i, j), d);
            for (j2 = j; j2 < A->c; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, i, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, i - 1, j2));
                fmpz_mul(fmpz_mat_entry(H, i - 1, j2), r1d,
                        fmpz_mat_entry(H, i - 1, j2));
                fmpz_submul(fmpz_mat_entry(H, i - 1, j2), r2d,
                        fmpz_mat_entry(H, i, j2));
                fmpz_set(fmpz_mat_entry(H, i, j2), b);
            }
        }
        fmpz_mat_swap_rows(H, NULL, A->r - 1, k);
        if (fmpz_sgn(fmpz_mat_entry(H, k, j)) < 0)
        {
            for (j2 = j; j2 < A->c; j2++)
            {
                fmpz_neg(fmpz_mat_entry(H, k, j2),
                        fmpz_mat_entry(H, k, j2));
            }
        }
        if (fmpz_is_zero(fmpz_mat_entry(H, k, j)))
        {
            k--;
            if (l > 0)
                l--;
        }
        else
        {
            /* reduce higher entries of column j with row k */
            for (i = k - 1; i >= 0; i--)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                        fmpz_mat_entry(H, k, j));
                for (j2 = j; j2 < A->c; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                            fmpz_mat_entry(H, k, j2));
                }
            }
        }
    }
    fmpz_clear(q);
    fmpz_clear(r2d);
    fmpz_clear(r1d);
    fmpz_clear(b);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(d);
}
