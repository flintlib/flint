/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void fmpz_mat_snf_kannan_bachem(fmpz_mat_t S, const fmpz_mat_t A)
{
    slong i, j, k, d, m, n;
    fmpz_t r1g, r2g, b, u, v, g;
    m = A->r;
    n = A->c;
    d = FLINT_MIN(m, n);
    fmpz_init(r1g);
    fmpz_init(r2g);
    fmpz_init(b);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(g);

    fmpz_mat_set(S, A);

    for (k = 0; k != d; k++)
    {
        int col_done;
        do
        {
            /* clear column */
            for (i = k + 1; i != m; i++)
            {
                /* reduce row i - 1 with row i */
                if (fmpz_is_zero(fmpz_mat_entry(S, i - 1, k)))
                    continue;
                if (fmpz_cmpabs(fmpz_mat_entry(S, i, k),
                            fmpz_mat_entry(S, i - 1, k)) == 0)
                {
                    if (fmpz_equal(fmpz_mat_entry(S, i, k),
                            fmpz_mat_entry(S, i - 1, k)))
                    {
                        for (j = k; j != n; j++)
                            fmpz_sub(fmpz_mat_entry(S, i - 1, j),
                                    fmpz_mat_entry(S, i - 1, j),
                                    fmpz_mat_entry(S, i, j));
                    }
                    else
                    {
                        for (j = k; j != n; j++)
                            fmpz_add(fmpz_mat_entry(S, i - 1, j),
                                    fmpz_mat_entry(S, i - 1, j),
                                    fmpz_mat_entry(S, i, j));
                    }
                    continue;
                }
                fmpz_xgcd(g, u, v, fmpz_mat_entry(S, i, k),
                        fmpz_mat_entry(S, i - 1, k));
                fmpz_divexact(r2g, fmpz_mat_entry(S, i - 1, k), g);
                fmpz_divexact(r1g, fmpz_mat_entry(S, i, k), g);
                for (j = k; j != n; j++)
                {
                    fmpz_mul(b, u, fmpz_mat_entry(S, i, j));
                    fmpz_addmul(b, v, fmpz_mat_entry(S, i - 1, j));
                    fmpz_mul(fmpz_mat_entry(S, i - 1, j), r1g,
                            fmpz_mat_entry(S, i - 1, j));
                    fmpz_submul(fmpz_mat_entry(S, i - 1, j), r2g,
                            fmpz_mat_entry(S, i, j));
                    fmpz_set(fmpz_mat_entry(S, i, j), b);
                }
            }
            fmpz_mat_swap_rows(S, NULL, m - 1, k);

            /* clear row */
            for (j = k + 1; j != n; j++)
            {
                /* reduce col j with col k */
                if (fmpz_is_zero(fmpz_mat_entry(S, k, j)))
                    continue;
                if (fmpz_cmpabs(fmpz_mat_entry(S, k, k),
                            fmpz_mat_entry(S, k, j)) == 0)
                {
                    if (fmpz_equal(fmpz_mat_entry(S, k, k),
                            fmpz_mat_entry(S, k, j)))
                    {
                        for (i = k; i != m; i++)
                            fmpz_sub(fmpz_mat_entry(S, i, j),
                                    fmpz_mat_entry(S, i, j),
                                    fmpz_mat_entry(S, i, k));
                    }
                    else
                    {
                        for (i = k; i != m; i++)
                            fmpz_add(fmpz_mat_entry(S, i, j),
                                    fmpz_mat_entry(S, i, j),
                                    fmpz_mat_entry(S, i, k));
                    }
                    continue;
                }
                fmpz_xgcd(g, u, v, fmpz_mat_entry(S, k, k),
                        fmpz_mat_entry(S, k, j));
                fmpz_divexact(r2g, fmpz_mat_entry(S, k, j), g);
                fmpz_divexact(r1g, fmpz_mat_entry(S, k, k), g);
                for (i = k; i != m; i++)
                {
                    fmpz_mul(b, u, fmpz_mat_entry(S, i, k));
                    fmpz_addmul(b, v, fmpz_mat_entry(S, i, j));
                    fmpz_mul(fmpz_mat_entry(S, i, j), r1g,
                            fmpz_mat_entry(S, i, j));
                    fmpz_submul(fmpz_mat_entry(S, i, j), r2g,
                            fmpz_mat_entry(S, i, k));
                    fmpz_set(fmpz_mat_entry(S, i, k), b);
                }
            }
            col_done = 1;
            for (i = 0; i != m; i++)
                col_done &= (i == k) || fmpz_is_zero(fmpz_mat_entry(S, i, k));
        }
        while (!col_done);

        if (fmpz_sgn(fmpz_mat_entry(S, k, k)) < 0)
            fmpz_neg(fmpz_mat_entry(S, k, k), fmpz_mat_entry(S, k, k));
    }

    fmpz_clear(r2g);
    fmpz_clear(r1g);
    fmpz_clear(b);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(g);

    fmpz_mat_snf_diagonal(S, S);
}
