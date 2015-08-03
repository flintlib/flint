/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "perm.h"

static void
add_columns(fmpz_mat_t H, const fmpz_mat_t B, const fmpz_mat_t H1)
{
    int neg;
    slong i, j, n, bits;
    fmpz_t den, tmp, one;
    fmpq_t num, alpha;
    fmpz_mat_t Bu, B1, cols, k;
    fmpq_mat_t H1_q, cols_q, x;
    flint_rand_t state;

    n = B->r;

    fmpz_mat_init(Bu, n, n);
    fmpz_mat_init(B1, n - 1, n);
    fmpz_mat_init(cols, n, B->c - n);
    fmpz_mat_init(k, n, 1);
    fmpq_mat_init(x, n, B->c - n);
    fmpq_mat_init(cols_q, n, B->c - n);
    fmpq_mat_init(H1_q, n, n);

    for (i = 0; i < n; i++)
        for (j = 0; j < cols->c; j++)
            fmpz_set(fmpz_mat_entry(cols, i, j), fmpz_mat_entry(B, i, n + j));
    for (i = 0; i < n - 1; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpz_set(fmpz_mat_entry(Bu, i, j), fmpz_mat_entry(B, i, j));
            fmpz_set(fmpz_mat_entry(B1, i, j), fmpz_mat_entry(B, i, j));
        }
    }

    /* find kernel basis vector */
    if (fmpz_mat_nullspace(k, B1) != 1)
    {
        flint_printf("Exception (fmpz_mat_hnf_pernet_stein). "
                     "Nullspace was not dimension one.\n");
        abort();
    }

    bits = fmpz_mat_max_bits(B1);
    if (bits < 0)
        bits = -bits;

    fmpz_mat_clear(B1);

    fmpz_init(tmp);

    /* set the last row of Bu to be random, such that Bu is nonsingular */
    flint_randinit(state);
    while (fmpz_is_zero(tmp))
    {
        _fmpz_vec_randtest(Bu->rows[n - 1], state, n, bits);

        fmpz_zero(tmp);
        for (j = 0; j < n; j++)
            fmpz_addmul(tmp, fmpz_mat_entry(Bu, n - 1, j),
                        fmpz_mat_entry(k, j, 0));
    }
    flint_randclear(state);
    fmpz_clear(tmp);

    /* solve Bu*x = cols */
    if (!fmpq_mat_solve_fmpz_mat(x, Bu, cols))
    {
        flint_printf("Exception (fmpz_mat_hnf_pernet_stein). "
                "Singular input matrix for solve.");
        abort();
    }

    /* fix final row */
    fmpq_init(num);
    fmpz_init(den);
    fmpq_init(alpha);
    fmpz_init(one);
    fmpz_one(one);

    /* compute denominator */
    for (i = 0; i < n; i++)
        fmpz_addmul(den, fmpz_mat_entry(B, n - 1, i), fmpz_mat_entry(k, i, 0));
    neg = (fmpz_sgn(den) < 0);
    if (neg)
        fmpz_neg(den, den);

    for (j = 0; j < B->c - H1->c; j++)
    {
        fmpq_zero(num);
        for (i = 0; i < n; i++)
        {
            _fmpq_addmul(fmpq_numref(num), fmpq_denref(num),
                         fmpz_mat_entry(B, n - 1, i), one,
                         fmpq_mat_entry_num(x, i, j),
                         fmpq_mat_entry_den(x, i, j));
        }
        _fmpq_sub(fmpq_numref(alpha), fmpq_denref(alpha),
                  fmpz_mat_entry(B, n - 1, n + j), one,
                  fmpq_numref(num), fmpq_denref(num));

        _fmpq_mul(fmpq_numref(alpha), fmpq_denref(alpha),
                  fmpq_numref(alpha), fmpq_denref(alpha), one, den);
        if (neg)
            fmpq_neg(alpha, alpha);

        /* x_i += alpha*k */
        for (i = 0; i < n; i++)
        {
            _fmpq_addmul(fmpq_mat_entry_num(x, i, j),
                         fmpq_mat_entry_den(x, i, j), fmpq_numref(alpha),
                         fmpq_denref(alpha), fmpz_mat_entry(k, i, 0), one);
        }
    }

    fmpq_clear(num);
    fmpz_clear(den);
    fmpz_clear(one);
    fmpq_clear(alpha);

    /* set cols = H1*x and place in position in H */
    fmpq_mat_set_fmpz_mat(H1_q, H1);
    fmpq_mat_mul(cols_q, H1_q, x);
    fmpq_mat_get_fmpz_mat(cols, cols_q);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H1, i, j));
        for (j = n; j < H->c; j++)
            fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(cols, i, j - n));
    }

    fmpq_mat_clear(H1_q);
    fmpq_mat_clear(x);
    fmpq_mat_clear(cols_q);
    fmpz_mat_clear(k);
    fmpz_mat_clear(cols);
    fmpz_mat_clear(Bu);
}

/* takes input matrix H with rows 0 to start_row - 1 in HNF to a HNF matrix */
static void
add_rows(fmpz_mat_t H, slong start_row)
{
    slong i, i2, j, j2, num_pivots, new_row, row;
    slong *pivots;
    fmpz_t b, d, u, v, r1d, r2d, q;

    fmpz_init(b);
    fmpz_init(d);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(r1d);
    fmpz_init(r2d);
    fmpz_init(q);
    num_pivots = 0;

    /* find the pivots of H */
    pivots = flint_malloc(H->c * sizeof(slong));
    for (i = j = 0; i < start_row; i++, j++)
    {
        for (; j < H->c && fmpz_is_zero(fmpz_mat_entry(H, i, j)); j++) ;
        if (j == H->c)
            break;
        pivots[i] = j;
        num_pivots = i + 1;
    }

    for (row = start_row; row < H->r; row++)
    {
        /* reduce row to be added with existing */
        for (i = 0; i < num_pivots; i++)
        {
            j = pivots[i];
            if (fmpz_is_zero(fmpz_mat_entry(H, row, j)))
                continue;
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, i, j),
                      fmpz_mat_entry(H, row, j));
            fmpz_divexact(r1d, fmpz_mat_entry(H, i, j), d);
            fmpz_divexact(r2d, fmpz_mat_entry(H, row, j), d);
            for (j2 = j; j2 < H->c; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, i, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, row, j2));
                fmpz_mul(fmpz_mat_entry(H, row, j2), r1d,
                         fmpz_mat_entry(H, row, j2));
                fmpz_submul(fmpz_mat_entry(H, row, j2), r2d,
                            fmpz_mat_entry(H, i, j2));
                fmpz_set(fmpz_mat_entry(H, i, j2), b);
            }
        }
        /* find first non-zero entry of the added row */
        for (j = 0; j < H->c && fmpz_is_zero(fmpz_mat_entry(H, row, j)); j++) ;
        new_row = row;
        if (j != H->c)          /* last row non-zero, move to correct position */
        {
            if (fmpz_sgn(fmpz_mat_entry(H, row, j)) < 0)
            {
                for (j2 = j; j2 < H->c; j2++)
                {
                    fmpz_neg(fmpz_mat_entry(H, row, j2),
                             fmpz_mat_entry(H, row, j2));
                }
            }
            do
            {
                if (new_row < row)
                    fmpz_mat_swap_rows(H, NULL, new_row, new_row + 1);
                new_row--;
                for (j2 = 0; j2 < H->c &&
                     fmpz_is_zero(fmpz_mat_entry(H, new_row, j2)); j2++) ;
            }
            while (j2 > j);
        }

        for (i = new_row, j = 0; i <= row; i++, j++)
        {
            for (; j < H->c && fmpz_is_zero(fmpz_mat_entry(H, i, j)); j++) ;
            if (j == H->c)
                break;
            pivots[i] = j;
            num_pivots = i + 1;
        }

        /* reduce above pivot entries */
        for (i = 0; i < num_pivots; i++)
        {
            for (i2 = 0; i2 < i; i2++)
            {
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i2, pivots[i]),
                            fmpz_mat_entry(H, i, pivots[i]));
                for (j2 = pivots[i]; j2 < H->c; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i2, j2), q,
                                fmpz_mat_entry(H, i, j2));
                }
            }
        }
    }

    fmpz_clear(q);
    fmpz_clear(r2d);
    fmpz_clear(r1d);
    fmpz_clear(v);
    fmpz_clear(u);
    fmpz_clear(d);
    fmpz_clear(b);
    flint_free(pivots);
}

static void
double_det(fmpz_t d1, fmpz_t d2, const fmpz_mat_t B, const fmpz_mat_t c,
           const fmpz_mat_t d)
{
    slong i, j, n;
    slong *P;
    mp_limb_t p, u1mod, u2mod, v1mod, v2mod;
    fmpz_t bound, prod, s1, s2, t, u1, u2, v1, v2;
    fmpz_mat_t dt, Bt;
    fmpq_t tmpq;
    fmpq_mat_t x;
    nmod_mat_t Btmod;

    n = B->c;

    fmpz_mat_init(dt, n, 1);
    fmpz_mat_init(Bt, n, n);
    fmpq_mat_init(x, n, 1);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - 1; j++)
            fmpz_set(fmpz_mat_entry(Bt, i, j), fmpz_mat_entry(B, j, i));
        fmpz_set(fmpz_mat_entry(Bt, i, n - 1), fmpz_mat_entry(c, 0, i));
    }

    /* solve B^Tx = d^T */
    fmpz_mat_transpose(dt, d);
    fmpq_mat_solve_fmpz_mat(x, Bt, dt);

    if (!fmpq_is_zero(fmpq_mat_entry(x, n - 1, 0)))
    {
        fmpz_init(bound);
        fmpz_init(prod);
        fmpz_init(t);
        fmpz_init(s1);
        fmpz_init(s2);
        fmpz_init(u1);
        fmpz_init(u2);
        fmpz_init(v1);
        fmpz_init(v2);

        /* compute lcm of denominators of vectors x and y */
        fmpq_init(tmpq);
        fmpz_one(u1);
        fmpz_one(u2);
        for (i = 0; i < n - 1; i++)
        {
            fmpz_lcm(u1, u1, fmpq_mat_entry_den(x, i, 0));
            fmpq_div(tmpq, fmpq_mat_entry(x, i, 0),
                     fmpq_mat_entry(x, n - 1, 0));
            fmpz_lcm(u2, u2, fmpq_denref(tmpq));
        }
        fmpz_lcm(u1, u1, fmpq_mat_entry_den(x, n - 1, 0));
        fmpq_inv(tmpq, fmpq_mat_entry(x, n - 1, 0));
        fmpz_lcm(u2, u2, fmpq_denref(tmpq));
        fmpq_clear(tmpq);

        /* compute Hadamard bounds */
        fmpz_one(bound);
        for (j = 0; j < n - 1; j++)
        {
            fmpz_zero(s1);
            for (i = 0; i < n; i++)
                fmpz_addmul(s1, fmpz_mat_entry(Bt, i, j),
                            fmpz_mat_entry(Bt, i, j));
            fmpz_sqrtrem(s1, t, s1);
            if (!fmpz_is_zero(t))
                fmpz_add_ui(s1, s1, UWORD(1));
            fmpz_mul(bound, bound, s1);
        }
        fmpz_zero(s1);
        fmpz_zero(s2);
        for (j = 0; j < n; j++)
        {
            fmpz_addmul(s1, fmpz_mat_entry(c, 0, j), fmpz_mat_entry(c, 0, j));
            fmpz_addmul(s2, fmpz_mat_entry(d, 0, j), fmpz_mat_entry(d, 0, j));
        }
        fmpz_sqrtrem(s1, t, s1);
        if (!fmpz_is_zero(t))
            fmpz_add_ui(s1, s1, UWORD(1));
        fmpz_sqrtrem(s2, t, s2);
        if (!fmpz_is_zero(t))
            fmpz_add_ui(s2, s2, UWORD(1));
        fmpz_mul(s1, s1, bound);
        fmpz_mul(s2, s2, bound);
        fmpz_cdiv_q(s1, s1, u1);
        fmpz_cdiv_q(s2, s2, u2);
        if (fmpz_cmp(s1, s2) > 0)
            fmpz_set(bound, s1);
        else
            fmpz_set(bound, s2);
        fmpz_mul_ui(bound, bound, UWORD(2));

        fmpz_one(prod);
        P = _perm_init(n);
        nmod_mat_init(Btmod, n, n, 2);
        p = UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS;
        /* compute determinants divided by u1 and u2 */
        while (fmpz_cmp(prod, bound) <= 0)
        {
            p = n_nextprime(p, 0);
            u1mod = fmpz_fdiv_ui(u1, p);
            u2mod = fmpz_fdiv_ui(u2, p);
            if (!(u1mod || u2mod))
                continue;
            _nmod_mat_set_mod(Btmod, p);
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n - 1; j++)
                    nmod_mat_entry(Btmod, i, j) =
                        fmpz_fdiv_ui(fmpz_mat_entry(B, j, i), p);
                nmod_mat_entry(Btmod, i, n - 1) =
                    fmpz_fdiv_ui(fmpz_mat_entry(c, 0, i), p);
            }
            nmod_mat_lu(P, Btmod, 0);
            v1mod = UWORD(1);
            for (i = 0; i < n; i++)
                v1mod = n_mulmod2_preinv(v1mod, nmod_mat_entry(Btmod, i, i), p,
                                         Btmod->mod.ninv);
            if (_perm_parity(P, n) == 1)
                v1mod = nmod_neg(v1mod, Btmod->mod);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n - 1; j++)
                    nmod_mat_entry(Btmod, i, j) =
                        fmpz_fdiv_ui(fmpz_mat_entry(B, j, i), p);
                nmod_mat_entry(Btmod, i, n - 1) =
                    fmpz_fdiv_ui(fmpz_mat_entry(d, 0, i), p);
            }
            nmod_mat_lu(P, Btmod, 0);
            v2mod = UWORD(1);
            for (i = 0; i < n; i++)
                v2mod = n_mulmod2_preinv(v2mod, nmod_mat_entry(Btmod, i, i), p,
                                         Btmod->mod.ninv);
            if (_perm_parity(P, n) == 1)
                v2mod = nmod_neg(v2mod, Btmod->mod);

            v1mod = n_mulmod2_preinv(v1mod, n_invmod(u1mod, p), p,
                                     Btmod->mod.ninv);
            v2mod = n_mulmod2_preinv(v2mod, n_invmod(u2mod, p), p,
                                     Btmod->mod.ninv);
            fmpz_CRT_ui(v1, v1, prod, v1mod, p, 1);
            fmpz_CRT_ui(v2, v2, prod, v2mod, p, 1);
            fmpz_mul_ui(prod, prod, p);
        }

        fmpz_mul(d1, u1, v1);
        fmpz_mul(d2, u2, v2);

        fmpz_clear(bound);
        fmpz_clear(prod);
        fmpz_clear(s1);
        fmpz_clear(s2);
        fmpz_clear(u1);
        fmpz_clear(u2);
        fmpz_clear(v1);
        fmpz_clear(v2);
        fmpz_clear(t);
        _perm_clear(P);
        nmod_mat_clear(Btmod);
    }
    else                        /* can't use the clever method above so naively compute both dets */
    {
        fmpz_mat_det(d1, Bt);
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(Bt, j, n - 1), fmpz_mat_entry(d, 0, j));
        fmpz_mat_det(d2, Bt);
    }

    fmpz_mat_clear(dt);
    fmpz_mat_clear(Bt);
    fmpq_mat_clear(x);
}

int
fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A, flint_rand_t state)
{
    slong i, j, m, n, p, r, *P, *pivots;
    fmpz_t d1, d2, g, s, t;
    fmpz_mat_t c, d, B, C, H1, H2, H3;
    nmod_mat_t Amod;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    /* find permutation so we can ensure first rows of H are nonsingular */
    P = _perm_init(m);
    pivots = _perm_init(n);
    
    p = n_randprime(state, NMOD_MAT_OPTIMAL_MODULUS_BITS, 1);
    nmod_mat_init(Amod, m, n, p);

    fmpz_mat_get_nmod_mat(Amod, A);
    r = _nmod_mat_rref(Amod, pivots, P);

    nmod_mat_clear(Amod);

    /* rank is not full, so silently fail */
    if (r < FLINT_MIN(m, n))
    {
        _perm_clear(P);
        _perm_clear(pivots);
        
        return 0;
    }

    /* if A has full column rank we might wish to use minors based hnf */
    if (r == n && n < 52)
    {
        slong b = fmpz_mat_max_bits(A), cutoff = 52;
        if (b < 0)
            b = -b;

        if (b <= 8)
            cutoff = 35;
        else if (b <= 32)
            cutoff = 44;
        else if (b <= 256)
            cutoff = 48;

        if (n < cutoff)
        {
            fmpz_mat_hnf_minors(H, A);
            _perm_clear(P);
            _perm_clear(pivots);
            return 1;
        }
    }

    fmpz_mat_init(c, 1, r - 1);
    fmpz_mat_init(d, 1, r - 1);
    fmpz_mat_init(B, r - 2, r - 1);
    fmpz_mat_init(C, r - 1, r - 1);

    for (i = 0; i < r - 2; i++)
    {
        for (j = 0; j < r - 1; j++)
        {
            fmpz_set(fmpz_mat_entry(B, i, j),
                     fmpz_mat_entry(A, P[i], pivots[j]));
            fmpz_set(fmpz_mat_entry(C, i, j),
                     fmpz_mat_entry(A, P[i], pivots[j]));
        }
        fmpz_set(fmpz_mat_entry(C, i, r - 1),
                 fmpz_mat_entry(A, P[i], pivots[r - 1]));
    }
    for (j = 0; j < r - 1; j++)
    {
        fmpz_set(fmpz_mat_entry(c, 0, j),
                 fmpz_mat_entry(A, P[r - 2], pivots[j]));
        fmpz_set(fmpz_mat_entry(d, 0, j),
                 fmpz_mat_entry(A, P[r - 1], pivots[j]));
    }

    fmpz_init(g);
    fmpz_init(s);
    fmpz_init(t);

    /* if rank is too low leave g = 0 so we don't try to decompose later */
    if (r > 2)
    {
        fmpz_init(d1);
        fmpz_init(d2);

        double_det(d1, d2, B, c, d);
        fmpz_xgcd(g, s, t, d1, d2);

        for (j = 0; j < r - 1; j++)
        {
            fmpz_mul(fmpz_mat_entry(C, r - 2, j), s,
                     fmpz_mat_entry(A, P[r - 2], pivots[j]));
            fmpz_addmul(fmpz_mat_entry(C, r - 2, j), t,
                        fmpz_mat_entry(A, P[r - 1], pivots[j]));
        }

        fmpz_clear(d2);
        fmpz_clear(d1);
    }

    if (!fmpz_is_zero(g))       /* chosen matrix invertible */
    {
        fmpz_mat_init(H1, r - 1, r - 1);

        if (COEFF_IS_MPZ(*g) && C->r > 3)   /* if g is too big, call another algorithm */
        {
            fmpz_mat_hnf(H1, C);
        }  else                    /* use modulo determinant algorithm to compute HNF of C */
            fmpz_mat_hnf_modular(H1, C, g);

        fmpz_mat_clear(B);
        fmpz_mat_init(B, r - 1, n);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < r - 2; i++)
                fmpz_set(fmpz_mat_entry(B, i, j),
                         fmpz_mat_entry(A, P[i], pivots[j]));
            fmpz_mul(fmpz_mat_entry(B, r - 2, j), s,
                     fmpz_mat_entry(A, P[r - 2], pivots[j]));
            fmpz_addmul(fmpz_mat_entry(B, r - 2, j), t,
                        fmpz_mat_entry(A, P[r - 1], pivots[j]));
        }

        fmpz_mat_init(H2, r - 1, n);
        fmpz_mat_init(H3, m + 1, n);

        add_columns(H2, B, H1);

        for (i = 0; i < r - 1; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(H3, i, pivots[j]),
                         fmpz_mat_entry(H2, i, j));

        for (i = 1; i <= m - r + 2; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(H3, H3->r - i, j),
                         fmpz_mat_entry(A, P[m - i], j));

        add_rows(H3, r - 1);

        /* fill H with HNF */
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H3, i, j));

        fmpz_mat_clear(H1);
        fmpz_mat_clear(H2);
        fmpz_mat_clear(H3);
    }
    else
    {
        if (r == n)             /* if A has full column rank we can use minors based hnf */
            fmpz_mat_hnf_minors(H, A);
        else
            fmpz_mat_hnf_classical(H, A);
    }

    _perm_clear(P);
    _perm_clear(pivots);
    fmpz_clear(t);
    fmpz_clear(s);
    fmpz_clear(g);
    fmpz_mat_clear(C);
    fmpz_mat_clear(B);
    fmpz_mat_clear(c);
    fmpz_mat_clear(d);

    return 1;
}
