/*
    Copyright (C) 2015 Tommy Hofmann 

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

static int
_fmpz_mat_pivot(fmpz_mat_t A, slong start_row, slong col)
{
    slong j;
    fmpz * u;

    if (fmpz_is_zero(fmpz_mat_entry(A, start_row, col)) == 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (!fmpz_is_zero(fmpz_mat_entry(A, j, col)))
        {
            u = A->rows[j];
            A->rows[j] = A->rows[start_row];
            A->rows[start_row] = u;

            return -1;
        }
    }
    return 0;
}

static void
_fmpz_ppio(fmpz_t ppi, fmpz_t ppo, fmpz_t a, fmpz_t b)
{
    fmpz_t c, n, g;

    fmpz_init(c);
    fmpz_init(n);
    fmpz_init(g);

    fmpz_gcd(c, a, b);
    fmpz_divexact(n, a, c);
    fmpz_gcd(g, c, n);

    while (!fmpz_is_one(g))
    {
        fmpz_mul(c, c, g);
        fmpz_divexact(n, n, g);
        fmpz_gcd(g, c, n);
    }
    fmpz_set(ppi, c);
    fmpz_set(ppo, n);

    fmpz_clear(c);
    fmpz_clear(n);
    fmpz_clear(g);
}

static void 
_fmpz_stab(fmpz_t t, const fmpz_t a, const fmpz_t b, const fmpz_t N)
{
    fmpz_t g, gg, s, aa, bb;

    fmpz_init(g);
    fmpz_init(gg);
    fmpz_init(s);
    fmpz_init(aa);
    fmpz_init(bb);

    fmpz_gcd(g, a, b);
    fmpz_gcd(gg, g, N);

    fmpz_divexact(bb, N, gg);

    fmpz_divexact(aa, a, gg);

    _fmpz_ppio(s, t, bb, aa);

    fmpz_clear(g);
    fmpz_clear(gg);
    fmpz_clear(s);
    fmpz_clear(aa);
    fmpz_clear(bb);
}

static void
_fmpz_unit(fmpz_t u, fmpz_t a, const fmpz_t N)
{
    fmpz_t g, s, t, l, d, k;

    fmpz_init(g);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_init(l);
    fmpz_init(d);
    fmpz_init(k);

    fmpz_xgcd(g, s, t, a, N);

    if (fmpz_is_one(g) == 1)
    {
        fmpz_set(u, s);
    }
    else
    {
        fmpz_divexact(l, N, g);
        _fmpz_stab(d, s, l, N);
        fmpz_mul(u, d, l);
        fmpz_add(u, u, s);
        fmpz_mod(u, u, N);
    }

    fmpz_clear(g);
    fmpz_clear(s);
    fmpz_clear(t);
    fmpz_clear(l);
    fmpz_clear(d);
    fmpz_clear(k);
}

/* test wether q*a = b mod N has a solution */
static int
_fmpz_is_divisible_mod(fmpz_t q, fmpz_t b, fmpz_t a, const fmpz_t N)
{
    fmpz_t g, e, t;

    fmpz_init(g);
    fmpz_init(e);
    fmpz_init(t);

    fmpz_xgcd(g, e, t, a, N);

    if (fmpz_divisible(b, g))
    {
        fmpz_divexact(g, b, g);
        fmpz_mul(g, e, g);
        fmpz_mod(q, g, N);

        fmpz_clear(g);
        fmpz_clear(e);
        fmpz_clear(t);
        return 1;
    }

    fmpz_clear(g);
    fmpz_clear(e);
    fmpz_clear(t);

    return 0;
}

void
fmpz_mat_strong_echelon_form_mod(fmpz_mat_t A, const fmpz_t mod)
{
    fmpz_t s, t, q, u, v, t1, t2, g;
    slong m, n, row, col, i, k, l;
    fmpz  ** r;
    fmpz * extra_row;

    if (fmpz_mat_is_empty(A))
        return;

    fmpz_init(s);
    fmpz_init(t);
    fmpz_init(q);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(g);

    n = A->r;
    m = A->c;
    r = A->rows;

    extra_row = _fmpz_vec_init(m);

    row = col = 0;

    for (row = 0; row < n; row++)
    {
        for (col = 0; col < m; col++)
        {
            fmpz_mod(fmpz_mat_entry(A, row, col), fmpz_mat_entry(A, row, col), mod);
        }
    }

    row = col = 0;

    while (row < n && col < m)
    {
        if (_fmpz_mat_pivot(A, row, col) == 0)
        {
            col++;
            continue;
        }
        for (i = row + 1; i < n; i++)
        {
            if (fmpz_is_zero(fmpz_mat_entry(A, i, col)))
            {
                continue;
            }
            if ( _fmpz_is_divisible_mod(s, fmpz_mat_entry(A, i, col), fmpz_mat_entry(A, row, col), mod))
            {
                for (k = col; k < m; k++)
                {
                    fmpz_set(t1, fmpz_mat_entry(A, i, k));
                    fmpz_submul(t1, s, fmpz_mat_entry(A, row, k));
                    fmpz_mod(t1, t1, mod);
                    fmpz_set(fmpz_mat_entry(A, i, k), t1);
                }
            }
            else
            {
                fmpz_xgcd(g, s, t, fmpz_mat_entry(A, row, col), fmpz_mat_entry(A, i, col));
                fmpz_divexact(u, fmpz_mat_entry(A, i, col), g);
                fmpz_neg(u, u);
                fmpz_divexact(v, fmpz_mat_entry(A, row, col), g);

                for (k = col; k < m; k++)
                {
                    fmpz_mul(t1, s, fmpz_mat_entry(A, row, k));
                    fmpz_addmul(t1, t, fmpz_mat_entry(A, i, k));
                    fmpz_mod(t1, t1, mod);
                    fmpz_mul(t2, u, fmpz_mat_entry(A, row, k));
                    fmpz_addmul(t2, v, fmpz_mat_entry(A, i, k));
                    fmpz_mod(t2, t2, mod);
                    fmpz_set(fmpz_mat_entry(A, row, k), t1);
                    fmpz_set(fmpz_mat_entry(A, i, k), t2);
                }
            }
        }
        for (i = row - 1; i >= 0; i--)
        {
            fmpz_mod(fmpz_mat_entry(A, i, col), fmpz_mat_entry(A, i, col), mod);
        }
        row++;
        col++;
    }

    for (col = 0; col < m; col++)
    {
        if (fmpz_is_zero(fmpz_mat_entry(A, col, col)) != 1)
        {
            _fmpz_unit(u, fmpz_mat_entry(A, col, col), mod);

            for (k = col ; k < m; k++)
            {
                fmpz_mul(fmpz_mat_entry(A, col, k), u, fmpz_mat_entry(A, col, k));
                fmpz_mod(fmpz_mat_entry(A, col, k), fmpz_mat_entry(A, col, k), mod);
            }

            for (row = 0; row < col ; row++)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(A, row, col), fmpz_mat_entry(A, col, col));

                for (l = row; l< m; l++)
                {
                    fmpz_submul(fmpz_mat_entry(A, row, l), q, fmpz_mat_entry(A, col, l));
                    fmpz_mod(fmpz_mat_entry(A, row, l), fmpz_mat_entry(A, row, l), mod);
                }
            }
            fmpz_gcd(g, mod, fmpz_mat_entry(A, col, col));

            if (fmpz_is_one(g) == 1)
            {
                continue;

            }
            fmpz_divexact(g, mod, g);
            _fmpz_vec_scalar_mul_fmpz(extra_row, r[col], m, g);
            _fmpz_vec_scalar_mod_fmpz(extra_row, extra_row, m, mod);
        }
        else
        {
          _fmpz_vec_set(extra_row, r[col], m);
        }

        for (row = col + 1; row < m; row++)
        {
            fmpz_xgcd(g, s, t, fmpz_mat_entry(A, row, row), extra_row + row);
            if (fmpz_is_zero(g))
            {
                continue;
            }
            fmpz_divexact(u, extra_row + row, g);
            fmpz_neg(u, u);
            fmpz_divexact(v, fmpz_mat_entry(A, row, row), g);

            for (k = row; k < m; k++)
            {
                fmpz_mul(t1, s, fmpz_mat_entry(A, row, k));
                fmpz_addmul(t1, t, extra_row + k);
                fmpz_mod(t1, t1, mod);
                fmpz_mul(t2, u, fmpz_mat_entry(A, row, k));
                fmpz_addmul(t2, v, extra_row + k);
                fmpz_mod(t2, t2, mod);
                fmpz_set(fmpz_mat_entry(A, row, k), t1);
                fmpz_set(extra_row + k, t2);
            }
        }
        if (fmpz_is_zero(fmpz_mat_entry(A, col, col)) != 1)
        {
            for (i = col - 1; i >= 0; i--)
            {
                fmpz_mod(fmpz_mat_entry(A, i, col), fmpz_mat_entry(A, i, col), mod);
            }
        }
    }

    for (i = 1; i < m; i++)
    {
        if (!fmpz_is_zero(fmpz_mat_entry(A, i, i)))
        {
            for (k = i - 1; k >= 0; k--)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(A, k, i), fmpz_mat_entry(A, i, i));
                for (l = i; l < m; l++)
                {
                    fmpz_submul(fmpz_mat_entry(A, k, l), fmpz_mat_entry(A, i, l), q);
                    fmpz_mod(fmpz_mat_entry(A, k, l), fmpz_mat_entry(A, k, l), mod);
                }
            }
        }
    }

    _fmpz_vec_clear(extra_row, m);
    fmpz_clear(s);
    fmpz_clear(t);
    fmpz_clear(q);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(g);
}

