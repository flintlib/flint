/*
    Copyright (C) 2024 Joseph Grantham

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
_fmpz_mat_snf_transform_diagonal(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V, const fmpz_mat_t A)
{
    fmpz_t t, x, y;
    slong i, j, k, n = FLINT_MIN(A->r, A->c);

    fmpz_init(t);
    fmpz_init(x);
    fmpz_init(y);

    fmpz_mat_set(S, A);
    for (i = 0; i < n; i++)
    {
       if (fmpz_sgn(fmpz_mat_entry(S, i, i)) < 0)
        {
            fmpz_neg(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i));
            for (k = 0; k < A->r; k++)
                fmpz_neg(fmpz_mat_entry(U, i, k), fmpz_mat_entry(U, i, k));
        }
    }
    for (j = n - 1; j >= 0; j--)
    {
        for (i = 0; i < j; i++)
        {
            if (fmpz_equal(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i + 1, i + 1)))
                continue;
            fmpz_xgcd(d, x, y, fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i + 1, i + 1));
            fmpz_divexact(fmpz_mat_entry(S, i + 1, i + 1), fmpz_mat_entry(S, i + 1, i + 1), d);

            /* record row operations */
            for (k = 0; k < A->r; k++)
            {
                fmpz_addmul(fmpz_mat_entry(U, i, k), y, fmpz_mat_entry(U, i + 1, k));
                fmpz_submul(fmpz_mat_entry(U, i + 1, k), fmpz_mat_entry(S, i + 1, i + 1),
                        fmpz_mat_entry(U, i, k));
                fmpz_neg(fmpz_mat_entry(U, i + 1, k), fmpz_mat_entry(U, i + 1, k));
            }

            fmpz_divexact(y, fmpz_mat_entry(S, i, i), d);

            /* record column operations */
            for (k = 0; k < A->c; k++)
            {
                fmpz_addmul(fmpz_mat_entry(V, k, i + 1), x, fmpz_mat_entry(V, k, i));
                fmpz_submul(fmpz_mat_entry(V, k, i), y, fmpz_mat_entry(V, k, i + 1));
                fmpz_swap(fmpz_mat_entry(V, k, i + 1), fmpz_mat_entry(V, k, i));
            }

            fmpz_mul(fmpz_mat_entry(S, i + 1, i + 1), fmpz_mat_entry(S, i + 1, i + 1), fmpz_mat_entry(S, i, i));
            fmpz_set(fmpz_mat_entry(S, i, i), d);
        }
    }
    fmpz_clear(t); 
    fmpz_clear(x);
    fmpz_clear(y);
}

void
_fmpz_mat_snf_transform_kannen_bachem(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V, const fmpz_mat_t A)
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
    fmpz_mat_one(U);
    fmpz_mat_one(V);


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
                        for (j = k; j != n; j++)
                            fmpz_sub(fmpz_mat_entry(U, i - 1, j),
                                    fmpz_mat_entry(U, i - 1, j),
                                    fmpz_mat_entry(U, i, j));
                    }
                    else
                    {
                        for (j = k; j != n; j++)
                            fmpz_add(fmpz_mat_entry(S, i - 1, j),
                                    fmpz_mat_entry(S, i - 1, j),
                                    fmpz_mat_entry(S, i, j));
                        for (j = k; j != n; j++)
                            fmpz_add(fmpz_mat_entry(U, i - 1, j),
                                    fmpz_mat_entry(U, i - 1, j),
                                    fmpz_mat_entry(U, i, j));
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
                for (j = k; j != n; j++)
                {
                    fmpz_mul(b, u, fmpz_mat_entry(U, i, j));
                    fmpz_addmul(b, v, fmpz_mat_entry(U, i - 1, j));
                    fmpz_mul(fmpz_mat_entry(U, i - 1, j), r1g,
                            fmpz_mat_entry(U, i - 1, j));
                    fmpz_submul(fmpz_mat_entry(U, i - 1, j), r2g,
                            fmpz_mat_entry(U, i, j));
                    fmpz_set(fmpz_mat_entry(U, i, j), b);
                }
            }
            fmpz_mat_swap_rows(S, NULL, m - 1, k);
            fmpz_mat_swap_rows(U, NULL, m - 1, k);

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
                        for (i = k; i != m; i++)
                            fmpz_sub(fmpz_mat_entry(V, i, j),
                                    fmpz_mat_entry(V, i, j),
                                    fmpz_mat_entry(V, i, k));
                    }
                    else
                    {
                        for (i = k; i != m; i++)
                            fmpz_add(fmpz_mat_entry(S, i, j),
                                    fmpz_mat_entry(S, i, j),
                                    fmpz_mat_entry(S, i, k));
                        for (i = k; i != m; i++)
                            fmpz_add(fmpz_mat_entry(V, i, j),
                                    fmpz_mat_entry(V, i, j),
                                    fmpz_mat_entry(V, i, k));
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
                for (i = k; i != m; i++)
                {
                    fmpz_mul(b, u, fmpz_mat_entry(V, i, k));
                    fmpz_addmul(b, v, fmpz_mat_entry(V, i, j));
                    fmpz_mul(fmpz_mat_entry(V, i, j), r1g,
                            fmpz_mat_entry(V, i, j));
                    fmpz_submul(fmpz_mat_entry(V, i, j), r2g,
                            fmpz_mat_entry(V, i, k));
                    fmpz_set(fmpz_mat_entry(V, i, k), b);
                }
            }
            col_done = 1;
            for (i = 0; i != m; i++)
                col_done &= (i == k) || fmpz_is_zero(fmpz_mat_entry(S, i, k));
        }
        while (!col_done);

        if (fmpz_sgn(fmpz_mat_entry(S, k, k)) < 0)
        {
            fmpz_neg(fmpz_mat_entry(S, k, k), fmpz_mat_entry(S, k, k));
            for (i = 0; i != m; i++)
                fmpz_neg(fmpz_mat_entry(U, k, i), fmpz_mat_entry(U, k, i));
        }
    }

    fmpz_clear(r2g);
    fmpz_clear(r1g);
    fmpz_clear(b);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(g);

    _fmpz_mat_snf_transform_diagonal(S, U, V, S);
}

void
fmpz_mat_snf_transform(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V, const fmpz_mat_t A) 
{
    _fmpz_mat_snf_transform_kannen_bachem(S, U, V, A);
}