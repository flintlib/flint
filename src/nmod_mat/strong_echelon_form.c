/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

static inline int
_nmod_mat_pivot(nmod_mat_t A, slong start_row, slong col)
{
    slong j;
    mp_ptr u;

    if (nmod_mat_entry(A, start_row, col) != 0)
        return 1;

    for (j = start_row + 1; j < A->r; j++)
    {
        if (nmod_mat_entry(A, j, col) != 0)
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
_n_ppio(mp_ptr ppi, mp_ptr ppo, mp_limb_t a, mp_limb_t b)
{
    mp_limb_t c, n, g;

    c = n_gcd(a, b);
    n = a/c;
    g = n_gcd(c, n);
    while( g != 1 )
    {
        c = c * g;
        n = n/g;
        g = n_gcd(c, n);
    }
    *ppi = c;
    *ppo = n;
}

static mp_limb_t
_n_stab(mp_limb_t a, mp_limb_t b, nmod_t N)
{
    mp_limb_t g, s, t;
    g = n_gcd(a, b);
    b = n_gcd(g, N.n);
    _n_ppio(&s, &t, N.n/b, a/b);
    return t;
}

static mp_limb_t
_n_unit(mp_limb_t a, nmod_t N)
{
    mp_limb_t g, s, l, d;

    g = n_gcdinv(&s, a, N.n);

    if (g == 1)
    {
        return s;
    }
    else
    {
        l = N.n/g;
        d = _n_stab(s, l, N);
        return nmod_add(s, nmod_mul(d, l, N), N);
    }
}

/* test whether q*a = b mod N has a solution */
static int
_n_is_divisible(mp_ptr q, mp_limb_t b, mp_limb_t a, nmod_t N)
{
    mp_limb_t e, g;
    g = n_gcdinv(&e, a, N.n);

    if (( b % g ) == 0)
    {
        *q = nmod_mul(e, b/g, N);
        return 1;
    }

    return 0;
}

void
nmod_mat_strong_echelon_form(nmod_mat_t A)
{
    mp_limb_t s, t, u, v, q, t1, t2, g;
    slong m, n, row, col, i, k, l;
    mp_limb_t **r;
    nmod_t mod;
    mp_ptr extra_row;

    if (nmod_mat_is_empty(A))
        return;

    n = A->r;
    m = A->c;
    r = A->rows;
    mod = A->mod;

    extra_row = _nmod_vec_init(m);

    row = col = 0;

    while (row < n && col < m)
    {
        if (_nmod_mat_pivot(A, row, col) == 0)
        {
            col++;
            continue;
        }
        for (i = row + 1; i < n; i++)
        {
            if (nmod_mat_entry(A, i, col) == 0)
            {
                continue;
            }

            if (_n_is_divisible(&s, nmod_mat_entry(A, i, col), nmod_mat_entry(A, row, col), mod))
            {
                 for (k = col; k < m; k++)
                 {
                     t1 = nmod_sub(nmod_mat_entry(A, i, k), nmod_mul(s, nmod_mat_entry(A, row, k), mod), mod);
                     nmod_mat_entry(A, i, k) = t1;
                 }
            }
            else
            {
                if (nmod_mat_entry(A, row, col) >= nmod_mat_entry(A, i, col))
                {
                    g = n_xgcd(&s, &t, nmod_mat_entry(A, row, col), nmod_mat_entry(A, i, col));
                }
                else
                {
                    g = n_xgcd(&t, &s, nmod_mat_entry(A, i, col), nmod_mat_entry(A, row, col));
                }

                /* now g = a*x - b*y
                 a,b < x < mod.n */
                t = nmod_neg(t, mod);
                u = (nmod_mat_entry(A, i, col))/g;
                u = nmod_neg(u, mod);
                v = (nmod_mat_entry(A, row, col))/g;
                /* now g = a*x + b*y and 0 = sv - tu = 1 modulo mod.n */

                for (k = col; k < m; k++)
                {
                    t1 = nmod_add(nmod_mul(s, nmod_mat_entry(A, row, k), mod), nmod_mul(t, nmod_mat_entry(A, i, k), mod), mod);
                    t2 = nmod_add(nmod_mul(u, nmod_mat_entry(A, row, k), mod), nmod_mul(v, nmod_mat_entry(A, i, k), mod), mod);
                    nmod_mat_entry(A, row, k) = t1;
                    nmod_mat_entry(A, i, k) = t2;
                }
            }
        }
        row++;
        col++;
    }

    for (col = 0; col < m; col++)
    {
        if (nmod_mat_entry(A, col, col) != 0)
        {
            u = _n_unit(nmod_mat_entry(A, col, col), mod);
            for (k = col; k < m; k++)
            {
                nmod_mat_entry(A, col, k) = nmod_mul(u, nmod_mat_entry(A, col, k), mod);
            }
            for (row = 0; row < col ; row++)
            {

                q = nmod_mat_entry(A, row, col)/nmod_mat_entry(A, col, col);

                for (l = row; l< m; l++)
                {
                    s = nmod_sub(nmod_mat_entry(A, row, l), nmod_mul(q, nmod_mat_entry(A, col, l), mod), mod);
                    nmod_mat_entry(A, row, l) = s;
                }
            }

            g = n_gcd(mod.n, nmod_mat_entry(A, col, col));
            if (g == 1)
            {
                continue;
            }
            g = mod.n/g;
            _nmod_vec_scalar_mul_nmod(extra_row, r[col], m, g, mod);
        }
        else
        {
          _nmod_vec_set(extra_row, r[col], m);
        }

        for (row = col + 1; row < m; row++)
        {
            if(nmod_mat_entry(A, row, row) >= extra_row[row])
            {
                g = n_xgcd(&s, &t, nmod_mat_entry(A, row, row), extra_row[row]);
            }
            else
            {
                g = n_xgcd(&t, &s, extra_row[row], nmod_mat_entry(A, row, row));
            }
            if (g == 0)
            {
                continue;
            }
            t = nmod_neg(t, mod);
            u = extra_row[row]/g;
            u = nmod_neg(u, mod);
            v = (nmod_mat_entry(A, row, row))/g;
            /* now g = a*x + b*y and 0 = sv - tu = 1 modulo mod.n */

            for (k = row; k < m; k++)
            {
                t1 = nmod_add(nmod_mul(s, nmod_mat_entry(A, row, k), mod), nmod_mul(t, extra_row[k], mod), mod);
                t2 = nmod_add(nmod_mul(u, nmod_mat_entry(A, row, k), mod), nmod_mul(v, extra_row[k], mod), mod);
                nmod_mat_entry(A, row, k) = t1;
                extra_row[k] = t2;

            }
        }
    }
    _nmod_vec_clear(extra_row);
}
