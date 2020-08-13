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
fmpz_mat_hnf_modular(fmpz_mat_t H, const fmpz_mat_t A, const fmpz_t D)
{
    slong j, i, k, m, n;
    fmpz_t R, R2, d, u, v, r1d, r2d, b, q;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_init_set(R, D);
    fmpz_init(R2);
    fmpz_init(u);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(d);
    fmpz_init(r1d);
    fmpz_init(r2d);
    fmpz_init(b);
    fmpz_init(q);
    fmpz_mat_set(H, A);

    for (k = 0; k != n; k++)
    {
        fmpz_fdiv_q_2exp(R2, R, 1);

        if (fmpz_is_zero(fmpz_mat_entry(H, k, k)))
            fmpz_set(fmpz_mat_entry(H, k, k), R);
        for (i = k + 1; i != m; i++)
        {
            /* reduce row i with row k mod R */
            if (fmpz_is_zero(fmpz_mat_entry(H, i, k)))
                continue;
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, k, k),
                      fmpz_mat_entry(H, i, k));
            fmpz_divexact(r1d, fmpz_mat_entry(H, k, k), d);
            fmpz_divexact(r2d, fmpz_mat_entry(H, i, k), d);
            for (j = k; j < n; j++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, k, j));
                fmpz_addmul(b, v, fmpz_mat_entry(H, i, j));
                fmpz_mul(fmpz_mat_entry(H, i, j), r1d,
                         fmpz_mat_entry(H, i, j));
                fmpz_submul(fmpz_mat_entry(H, i, j), r2d,
                            fmpz_mat_entry(H, k, j));
                fmpz_mod(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H, i, j), R);
                if (fmpz_cmp(fmpz_mat_entry(H, i, j), R2) > 0)
                    fmpz_sub(fmpz_mat_entry(H, i, j),
                             fmpz_mat_entry(H, i, j), R);
                fmpz_mod(fmpz_mat_entry(H, k, j), b, R);
                if (fmpz_cmp(fmpz_mat_entry(H, k, j), R2) > 0)
                    fmpz_sub(fmpz_mat_entry(H, k, j),
                             fmpz_mat_entry(H, k, j), R);
            }
        }

        fmpz_xgcd(d, u, v, fmpz_mat_entry(H, k, k), R);
        for (j = k; j < n; j++)
        {
            fmpz_mul(fmpz_mat_entry(H, k, j), u, fmpz_mat_entry(H, k, j));
            fmpz_mod(fmpz_mat_entry(H, k, j), fmpz_mat_entry(H, k, j), R);
        }
        if (fmpz_is_zero(fmpz_mat_entry(H, k, k)))
            fmpz_set(fmpz_mat_entry(H, k, k), R);

        /* reduce higher entries of column k with row k */
        for (i = k - 1; i >= 0; i--)
        {
            fmpz_fdiv_q(q, fmpz_mat_entry(H, i, k), fmpz_mat_entry(H, k, k));
            for (j = k; j < n; j++)
            {
                fmpz_submul(fmpz_mat_entry(H, i, j), q,
                            fmpz_mat_entry(H, k, j));
            }
        }
        fmpz_divexact(R, R, d);
    }

    fmpz_clear(b);
    fmpz_clear(r2d);
    fmpz_clear(r1d);
    fmpz_clear(q);
    fmpz_clear(d);
    fmpz_clear(v);
    fmpz_clear(u);
    fmpz_clear(R2);
    fmpz_clear(R);
}
