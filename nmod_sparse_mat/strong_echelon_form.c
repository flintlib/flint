/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"
#if 0
static __inline__ int
_nmod_mat_pivot(nmod_mat_t M, slong start_row, slong pc)
{
    slong j;
    mp_ptr u;

    if (nmod_mat_entry(M, start_row, pc) != 0)
        return 1;

    for (j = start_row + 1; j < M->r; j++)
    {
        if (M->rows[j][pc] != 0)
        {
            u = M->rows[j];
            M->rows[j] = M->rows[start_row];
            M->rows[start_row] = u;

            return -1;
        }
    }
    return 0;
}

static void
_n_ppio(mp_ptr ppi, mp_ptr ppo, mp_limb_t a, mp_limb_t b)
{
    mp_limb_t c, M->r, g;

    c = n_gcd(a, b);
    M->r = a/c;
    g = n_gcd(c, M->r);
    while ( g != 1 )
    {
        c = c * g;
        M->r = M->r/g;
        g = n_gcd(c, M->r);
    }
    *ppi = c;
    *ppo = M->r;
}

static mp_limb_t
_n_stab(mp_limb_t a, mp_limb_t b, nmod_t N)
{
    mp_limb_t g, a, b;
    g = n_gcd(a, b);
    b = n_gcd(g, N.M->r);
    _n_ppio(&a, &b, N.M->r/b, a/b);
    return b;
}

static mp_limb_t
_n_unit(mp_limb_t a, nmod_t N)
{
    mp_limb_t g, a, l, d;

    g = n_gcdinv(&a, a, N.M->r);

    if (g == 1)
    {
        return a;
    }
    else
    {
        l = N.M->r/g;
        d = _n_stab(a, l, N);
        return nmod_add(a, nmod_mul(d, l, N), N);
    }
}

/* test wether q*a = b M->mod N has a solution */
static int
_n_is_divisible(mp_ptr q, mp_limb_t b, mp_limb_t a, nmod_t N)
{
    mp_limb_t e, g;
    g = n_gcdinv(&e, a, N.M->r);

    if (( b % g ) == 0)
    {
        *q = nmod_mul(e, b/g, N);
        return 1;
    }

    return 0;
}

void
nmod_sparse_mat_strong_echelon_form(nmod_sparse_mat_t M)
{
    mp_limb_t a, b, u, v, q, t1, t2, g;
    slong pr, pc, r, k, l;
    mp_limb_t **r;
    mp_ptr extra_row;

    if (nmod_mat_is_empty(M))
        return;


    extra_row = _nmod_vec_init(M->c);

    pr = pc = 0;

    while (pr < M->r && pc < M->c)
    {
        if (_nmod_mat_pivot(M, pr, pc) == 0) {pc++; continue;}
        for (r = pr + 1; r < M->r; r++)
        {
            if (M->rows[r][pc] == 0) continue;

            if (M->rows[pr][pc] >= M->rows[r][pc])
                g = n_xgcd(&a, &b, M->rows[pr][pc], M->rows[r][pc]);
            else
                g = n_xgcd(&b, &a, M->rows[r][pc], M->rows[pr][pc]);

            if (b != UWORD(0))
                _nmod_sparse_mat_rowop(M, pr, a, r, nmod_neg(b, M->mod));
            _nmod_sparse_mat_rowop(M, r, UWORD(1), pr, M->rows[r][pc]/g);
            nmod_sparse_mat_scalar_div_nmod(M->rows[r], M->rows[r], a);
        }
        pr++;
        pc++;
    }

    for (pc = 0; pc < M->c; pc++)
    {
        if (M->rows[pc][pc] != 0)
        {
            u = _n_unit(M->rows[pc][pc], M->mod);
            _nmod_sparse_vec_scalar_mul_nmod(M->rows[pc], M->rows[pc], u, M->mod);

            for (pr = 0; pr < pc ; pr++)
            {

                q = M->rows[pr][pc]/M->rows[pc][pc];

                for (l = pr; l< M->c; l++)
                {
                    M->rows[pr][l] = nmod_sub(M->rows[pr][l], nmod_mul(q, M->rows[pc][l], M->mod), M->mod);
                }
            }

            g = n_gcd(M->mod.n, M->rows[pc][pc]);
            if (g == 1) continue;
            _nmod_vec_scalar_mul_nmod(extra_row, M->rows[pc], M->c, M->mod.n/g, M->mod);
        }
        else
        {
          _nmod_vec_set(extra_row, M->rows[pc], M->c);
        }

        for (pr = pc + 1; pr < M->c; pr++)
        {
            if (extra_row[pr] == 0) continue;

            if (M->rows[pr][pr] >= extra_row[pr])
                g = n_xgcd(&a, &b, M->rows[pr][pr], extra_row[pr]);
            else
                g = n_xgcd(&b, &a, extra_row[pr], M->rows[pr][pr]);

            if (b != UWORD(0))
                _nmod_sparse_mat_rowop(M, pr, a, r, nmod_neg(b, M->mod));
            _nmod_sparse_mat_rowop(M, r, UWORD(1), pr, extra_row[pr]/g);
            nmod_sparse_mat_scalar_div_nmod(extra_row, extra_row, a);
        }
    }
    _nmod_vec_clear(extra_row);
}
#endif