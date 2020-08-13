/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

static void _eliminate_col(fmpz_mat_t S, slong i, const fmpz_t mod)
{
    slong j, k, m, n;
    fmpz * t;
    fmpz_t b, g, u, v, r1g, r2g;

    m = S->r;
    n = S->c;

    if (i == m - 1)
    {
        fmpz_gcd(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i), mod);
        return;
    }

    fmpz_init(g);
    fmpz_init(u);
    fmpz_init(b);
    fmpz_init(r1g);
    fmpz_init(r2g);

    if (!fmpz_is_zero(fmpz_mat_entry(S, i, i)))
    {
        fmpz_init(v);

        fmpz_xgcd(g, u, v, fmpz_mat_entry(S, i + 1, i),
                fmpz_mat_entry(S, i, i));
        fmpz_divexact(r1g, fmpz_mat_entry(S, i + 1, i), g);
        fmpz_divexact(r2g, fmpz_mat_entry(S, i, i), g);
        for (j = i; j < n; j++)
        {
            fmpz_mul(b, u, fmpz_mat_entry(S, i + 1, j));
            fmpz_addmul(b, v, fmpz_mat_entry(S, i, j));
            fmpz_mul(fmpz_mat_entry(S, i, j), r1g,
                    fmpz_mat_entry(S, i, j));
            fmpz_submul(fmpz_mat_entry(S, i, j), r2g,
                    fmpz_mat_entry(S, i + 1, j));
            fmpz_set(fmpz_mat_entry(S, i + 1, j), b);
        }

        fmpz_clear(v);
    }

    /* compute extended gcd of entries in column i */
    t = _fmpz_vec_init(m - i - 1);

    fmpz_set(g, fmpz_mat_entry(S, i + 1, i));
    fmpz_one(t);
    for (j = 2; j < m - i; j++)
    {
        fmpz_xgcd(g, u, t + j - 1, g, fmpz_mat_entry(S, i + j, i));
        for (k = 0; k < j - 1; k++)
            fmpz_mul(t + k, t + k, u);
    }

    /* set row i to have gcd in col i */
    for (k = i + 1; k < m; k++)
    {
        fmpz_mod(t + k - i - 1, t + k - i - 1, mod);
        for (j = i; j < n; j++)
            fmpz_addmul(fmpz_mat_entry(S, i, j), t + k - i - 1,
                    fmpz_mat_entry(S, k, j));
    }

    _fmpz_vec_clear(t, m - i - 1);

    /* reduce each row k with row i */
    if (!fmpz_is_zero(g)) /* if g = 0 then don't need to reduce */
    {
        for (k = i + 1; k < m; k++)
        {
            fmpz_divexact(r1g, fmpz_mat_entry(S, k, i), g);
            fmpz_neg(r1g, r1g);
            for (j = i; j < n; j++)
                fmpz_addmul(fmpz_mat_entry(S, k, j), r1g,
                        fmpz_mat_entry(S, i, j));
        }
        for (k = i + 1; k < m; k++)
            fmpz_mod(fmpz_mat_entry(S, k, i), fmpz_mat_entry(S, k, i), mod);
    }
    for (j = i; j < m; j++)
        for (k = i + 1; k < n; k++)
            fmpz_fdiv_r(fmpz_mat_entry(S, j, k), fmpz_mat_entry(S, j, k), mod);
    fmpz_gcd(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i), mod);

    fmpz_clear(b);
    fmpz_clear(g);
    fmpz_clear(u);
    fmpz_clear(r1g);
    fmpz_clear(r2g);
}

static void _eliminate_row(fmpz_mat_t S, slong i, const fmpz_t mod)
{
    slong j, k, m, n;
    fmpz * t;
    fmpz_t b, g, u, v, r1g, r2g, halfmod;

    m = S->r;
    n = S->c;

    if (i == n - 1)
    {
        fmpz_gcd(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i), mod);
        return;
    }

    fmpz_init(g);
    fmpz_init(u);
    fmpz_init(b);
    fmpz_init(r1g);
    fmpz_init(r2g);
    fmpz_init(halfmod);
    fmpz_fdiv_q_2exp(halfmod, mod, 1);

    if (!fmpz_is_zero(fmpz_mat_entry(S, i, i)))
    {
        fmpz_init(v);

        fmpz_xgcd(g, u, v, fmpz_mat_entry(S, i, i + 1),
                fmpz_mat_entry(S, i, i));
        fmpz_divexact(r1g, fmpz_mat_entry(S, i, i + 1), g);
        fmpz_divexact(r2g, fmpz_mat_entry(S, i, i), g);
        for (j = i; j < m; j++)
        {
            fmpz_mul(b, u, fmpz_mat_entry(S, j, i + 1));
            fmpz_addmul(b, v, fmpz_mat_entry(S, j, i));
            fmpz_mul(fmpz_mat_entry(S, j, i), r1g,
                    fmpz_mat_entry(S, j, i));
            fmpz_submul(fmpz_mat_entry(S, j, i), r2g,
                    fmpz_mat_entry(S, j, i + 1));
            fmpz_set(fmpz_mat_entry(S, j, i + 1), b);
        }

        fmpz_clear(v);
    }

    /* compute extended gcd of entries in row i */
    t = _fmpz_vec_init(n - i - 1);

    fmpz_set(g, fmpz_mat_entry(S, i, i + 1));
    fmpz_one(t);
    for (j = 2; j < n - i; j++)
    {
        fmpz_xgcd(g, u, t + j - 1, g, fmpz_mat_entry(S, i, i + j));
        for (k = 0; k < j - 1; k++)
            fmpz_mul(t + k, t + k, u);
    }

    /* reduce col i to have gcd in row i */
    for (k = i + 1; k < n; k++)
    {
        fmpz_mod(t + k - i - 1, t + k - i - 1, mod);
        for (j = i; j < m; j++)
            fmpz_addmul(fmpz_mat_entry(S, j, i), t + k - i - 1,
                    fmpz_mat_entry(S, j, k));
    }

    _fmpz_vec_clear(t, n - i - 1);

    /* reduce each col k with col i */
    if (!fmpz_is_zero(g)) /* if g = 0 then don't need to reduce */
    {
        for (k = i + 1; k < n; k++)
        {
            fmpz_divexact(r1g, fmpz_mat_entry(S, i, k), g);
            fmpz_neg(r1g, r1g);
            for (j = i; j < m; j++)
                fmpz_addmul(fmpz_mat_entry(S, j, k), r1g,
                        fmpz_mat_entry(S, j, i));
        }
    }
    for (j = i + 1; j < m; j++)
        for (k = i; k < n; k++)
            fmpz_fdiv_r(fmpz_mat_entry(S, j, k), fmpz_mat_entry(S, j, k), mod);
    fmpz_gcd(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i), mod);

    fmpz_clear(b);
    fmpz_clear(g);
    fmpz_clear(u);
    fmpz_clear(r1g);
    fmpz_clear(r2g);
    fmpz_clear(halfmod);
}

void fmpz_mat_snf_iliopoulos(fmpz_mat_t S, const fmpz_mat_t A, const fmpz_t mod)
{
    slong i, k, n;
    int done;

    n = FLINT_MIN(A->c, A->r);

    fmpz_mat_set(S, A);

    for (i = 0; i < A->r; i++)
        for (k = 0; k < A->c; k++)
            fmpz_mod(fmpz_mat_entry(S, i, k), fmpz_mat_entry(S, i, k), mod);

    for (k = 0; k != n; k++)
    {
        do
        {
            _eliminate_row(S, k, mod);
            _eliminate_col(S, k, mod);
            done = 1;
            if (fmpz_is_zero(fmpz_mat_entry(S, k, k)))
            {
                for (i = k + 1; i < A->c && done; i++)
                    done = fmpz_is_zero(fmpz_mat_entry(S, k, i));
            }
            else
            {
                for (i = k + 1; i < A->c && done; i++)
                    done = fmpz_divisible(fmpz_mat_entry(S, k, i),
                            fmpz_mat_entry(S, k, k));
            }
        }
        while (!done);
        for (i = k + 1; i < A->c; i++)
            fmpz_zero(fmpz_mat_entry(S, k, i));
    }

    fmpz_mat_snf_diagonal(S, S);
}
