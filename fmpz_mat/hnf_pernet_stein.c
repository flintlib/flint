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

void add_column(fmpz_mat_t H, const fmpz_mat_t B, const fmpz_mat_t H1)
{
    slong i, j, n, bits;
    fmpz_t M, den, one;
    fmpq_t num, alpha;
    fmpz_mat_t Bu, B1, col, x, k;
    fmpq_mat_t H1_q, col_q, x_q;
    flint_rand_t state;

    n = B->r;

    fmpz_init(M);
    fmpq_init(alpha);
    fmpz_init(one);
    fmpq_init(num);
    fmpz_init(den);
    fmpz_one(one);
    flint_randinit(state);
    fmpz_mat_init(Bu, n, n);
    fmpz_mat_init(B1, n - 1, n);
    fmpz_mat_init(col, n, 1);
    fmpz_mat_init(x, n, 1);
    fmpz_mat_init(k, n, 1);
    fmpq_mat_init(x_q, n, 1);
    fmpq_mat_init(col_q, n, 1);
    fmpq_mat_init(H1_q, H1->r, H1->c);

    for (i = 0; i < n; i++)
        fmpz_set(fmpz_mat_entry(col, i, 0), fmpz_mat_entry(B, i, n));
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
    fmpq_mat_set_fmpz_mat(H1_q, H1);
    while (1)
    {
        fmpz_zero(M);
        _fmpz_vec_randtest(Bu->rows[n - 1], state, Bu->c, bits);
        for (j = 0; j < Bu->c; j++)
            fmpz_addmul(M, fmpz_mat_entry(Bu, Bu->r - 1, j),
                    fmpz_mat_entry(k, j, 0));
        /*if (fmpz_mat_rank(Bu) != Bu->r) continue;*/
        if (fmpz_is_zero(M)) continue;
        /* TODO use good dixon code here */
        fmpz_mat_solve_dixon(x, M, Bu, col);
        fmpq_mat_set_fmpz_mat_mod_fmpz(x_q, x, M);
        fmpq_zero(num);
        fmpz_zero(den);
        for (i = 0; i < B->c - 1; i++)
        {
            _fmpq_addmul(fmpq_numref(num), fmpq_denref(num),
                    fmpz_mat_entry(B, n - 1, i), one,
                    fmpq_mat_entry_num(x_q, i, 0),
                    fmpq_mat_entry_den(x_q, i, 0));
            fmpz_addmul(den, fmpz_mat_entry(B, n - 1, i),
                    fmpz_mat_entry(k, i, 0));
        }
        if (0)
        {
            flint_printf("Exception (fmpz_mat_hnf_pernet_stein). "
                    "Thing does not divide thing.\n");
            fmpq_print(num); flint_printf("\n");
            fmpz_print(den);
            abort();
        }
        _fmpq_sub(fmpq_numref(alpha), fmpq_denref(alpha),
                fmpz_mat_entry(B, n - 1, n), one,
                fmpq_numref(num), fmpq_denref(num));
        if (fmpz_sgn(den) < 0)
        {
            fmpz_abs(den, den);
            _fmpq_mul(fmpq_numref(alpha), fmpq_denref(alpha),
                    fmpq_numref(alpha), fmpq_denref(alpha), one, den);
            fmpq_neg(alpha, alpha);
        }
        else
        {
            _fmpq_mul(fmpq_numref(alpha), fmpq_denref(alpha),
                    fmpq_numref(alpha), fmpq_denref(alpha), one, den);
        }
        /* x += alpha*k */
        for (i = 0; i < n; i++)
        {
            _fmpq_addmul(fmpq_mat_entry_num(x_q, i, 0),
                    fmpq_mat_entry_den(x_q, i, 0), fmpq_numref(alpha),
                    fmpq_denref(alpha), fmpz_mat_entry(k, i, 0), one);
        }
        fmpq_mat_mul(col_q, H1_q, x_q);
        fmpq_mat_get_fmpz_mat(col, col_q);
        for (i = 0; i < H1->r; i++)
        {
            for (j = 0; j < H1->c; j++)
                fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H1, i, j));
            fmpz_set(fmpz_mat_entry(H, i, H1->c), fmpz_mat_entry(col, i, 0));
        }
        break;
    }

    flint_randclear(state);
    fmpq_mat_clear(H1_q);
    fmpq_mat_clear(x_q);
    fmpq_mat_clear(col_q);
    fmpz_mat_clear(x);
    fmpz_mat_clear(k);
    fmpz_mat_clear(col);
    fmpz_mat_clear(Bu);
    fmpz_mat_clear(B1);
    fmpq_clear(num);
    fmpz_clear(den);
    fmpz_clear(one);
    fmpz_clear(M);
    fmpq_clear(alpha);
}

void add_row(fmpz_mat_t H, const fmpz_mat_t A, const fmpz * row)
{
    slong i, i2, j, j2, num_pivots, new_row;
    slong * pivots;
    fmpz_t b, d, u, v, r1d, r2d, q;

    fmpz_init(b);
    fmpz_init(d);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(r1d);
    fmpz_init(r2d);
    fmpz_init(q);
    num_pivots = 0;

    for (j = 0; j < A->c; j++)
    {
        for (i = 0; i < A->r; i++)
            fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(A, i, j));
        fmpz_set(fmpz_mat_entry(H, A->r, j), &row[j]);
    }

    /* find the pivots of A */
    for (i = j = 0; i < A->r; i++, j++)
    {
        for (; j < A->c && fmpz_is_zero(fmpz_mat_entry(A, i, j)); j++);
        if (j == A->c)
            break;
        num_pivots = i + 1;
    }
    pivots = flint_malloc(num_pivots * sizeof(slong));
    for (i = j = 0; i < A->r; i++, j++)
    {
        for (; j < A->c && fmpz_is_zero(fmpz_mat_entry(A, i, j)); j++);
        if (j == A->c)
            break;
        pivots[i] = j;
    }

    /* reduce row to be added with existing */
    for (i = 0; i < num_pivots; i++)
    {
        j = pivots[i];
        if (fmpz_is_zero(fmpz_mat_entry(H, A->r, j)))
            continue;
        fmpz_xgcd(d, u, v, fmpz_mat_entry(H, i, j), fmpz_mat_entry(H, A->r, j));
        fmpz_divexact(r1d, fmpz_mat_entry(H, i, j), d);
        fmpz_divexact(r2d, fmpz_mat_entry(H, A->r, j), d);
        for (j2 = j; j2 < A->c; j2++)
        {
            fmpz_mul(b, u, fmpz_mat_entry(H, i, j2));
            fmpz_addmul(b, v, fmpz_mat_entry(H, A->r, j2));
            fmpz_mul(fmpz_mat_entry(H, A->r, j2), r1d,
                    fmpz_mat_entry(H, A->r, j2));
            fmpz_submul(fmpz_mat_entry(H, A->r, j2), r2d,
                    fmpz_mat_entry(H, i, j2));
            fmpz_set(fmpz_mat_entry(H, i, j2), b);
        }
    }
    /* find first non-zero entry of the added row */
    for (j = 0; j < A->c && fmpz_is_zero(fmpz_mat_entry(H, A->r, j)); j++);
    new_row = A->r;
    if (j != A->c) /* last row non-zero, move to correct position */
    {
        if (fmpz_sgn(fmpz_mat_entry(H, A->r, j)) < 0)
        {
            for (j2 = j; j2 < A->c; j2++)
            {
                fmpz_neg(fmpz_mat_entry(H, A->r, j2),
                        fmpz_mat_entry(H, A->r, j2));
            }
        }
        do
        {
            if (new_row < A->r)
                fmpz_mat_swap_rows(H, NULL, new_row, new_row + 1);
            new_row--;
            for (j2 = 0; j2 < A->c &&
                    fmpz_is_zero(fmpz_mat_entry(H, new_row, j2)); j2++);
        }
        while (j2 > j);
    }

    for (i = j = num_pivots = 0; i < H->r; i++, j++)
    {
        for (; j < H->c && fmpz_is_zero(fmpz_mat_entry(H, i, j)); j++);
        if (j == H->c)
            break;
        num_pivots = i + 1;
    }
    pivots = flint_realloc(pivots, num_pivots * sizeof(slong));
    for (i = new_row, j = 0; i < H->r; i++, j++)
    {
        for (; j < H->c && fmpz_is_zero(fmpz_mat_entry(H, i, j)); j++);
        if (j == H->c)
            break;
        pivots[i] = j;
    }

    for (i = 0; i < num_pivots; i++)
    {
        for (i2 = 0; i2 < i; i2++)
        {
            fmpz_fdiv_q(q, fmpz_mat_entry(H, i2, pivots[i]),
                    fmpz_mat_entry(H, i, pivots[i]));
            for (j2 = pivots[i]; j2 < A->c; j2++)
            {
                fmpz_submul(fmpz_mat_entry(H, i2, j2), q,
                        fmpz_mat_entry(H, i, j2));
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

void double_det(fmpz_t d1, fmpz_t d2, const fmpz_mat_t B, const fmpz_mat_t c,
        const fmpz_mat_t d)
{
    slong i, j, n;
    slong * P;
    mp_limb_t q, p, u1mod, u2mod, v1mod, v2mod;
    fmpz_t bound, M, prod, s1, s2, t, u1, u2, v1, v2;
    fmpz_mat_t dt, Bt, x;
    fmpq_t tmp;
    fmpq_mat_t x_q;
    nmod_mat_t Bmod, Btmod, cmod, dmod;

    n = B->c;

    fmpz_mat_init(dt, n, 1);
    fmpz_mat_init(Bt, n, n);
    fmpz_mat_init(x, n, 1);
    fmpq_mat_init(x_q, n, 1);
    fmpz_init(M);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n - 1; j++)
            fmpz_set(fmpz_mat_entry(Bt, i, j), fmpz_mat_entry(B, j, i));
        fmpz_set(fmpz_mat_entry(Bt, i, n - 1), fmpz_mat_entry(c, 0, i));
    }
    fmpz_mat_transpose(dt, d);
    if (fmpz_mat_solve_dixon(x, M, Bt, dt))
        fmpq_mat_set_fmpz_mat_mod_fmpz(x_q, x, M);
    else
        fmpq_zero(fmpq_mat_entry(x_q, n - 1, 0));
    if (!fmpq_is_zero(fmpq_mat_entry(x_q, n - 1, 0)))
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
        fmpq_init(tmp);
        fmpz_one(u1);
        fmpz_one(u2);
        for (i = 0; i < n - 1; i++)
        {
            fmpz_lcm(u1, u1, fmpq_mat_entry_den(x_q, i, 0));
            fmpq_div(tmp, fmpq_mat_entry(x_q, i, 0),
                    fmpq_mat_entry(x_q, n - 1, 0));
            fmpz_lcm(u2, u2, fmpq_denref(tmp));
        }
        fmpz_lcm(u1, u1, fmpq_mat_entry_den(x_q, n - 1, 0));
        fmpq_inv(tmp, fmpq_mat_entry(x_q, n - 1, 0));
        fmpz_lcm(u2, u2, fmpq_denref(tmp));
        fmpq_clear(tmp);

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
        P = flint_malloc(n * sizeof(slong));
        nmod_mat_init(Btmod, n, n, 2);
        /*nmod_mat_init(Bmod, n, n - 1, 2);
        nmod_mat_init(cmod, n, 1, 2);
        nmod_mat_init(dmod, n, 1, 2);*/
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
            /*_nmod_mat_set_mod(Bmod, p);
            _nmod_mat_set_mod(cmod, p);
            _nmod_mat_set_mod(dmod, p);
            for (i = 0; i < n; i++)
            {
                nmod_mat_entry(cmod, i, 0) = fmpz_fdiv_ui(fmpz_mat_entry(c, 0, i), p);
                nmod_mat_entry(dmod, i, 0) = fmpz_fdiv_ui(fmpz_mat_entry(d, 0, i), p);
                for (j = 0; j < n - 1; j++)
                    nmod_mat_entry(Bmod, i, j) = fmpz_fdiv_ui(fmpz_mat_entry(B, j, i), p);
            }
            nmod_mat_lu(P, Bmod, 0);*/
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
            /*nmod_mat_solve_tril(cmod, Bmod, cmod, 1);
            nmod_mat_solve_tril(dmod, Bmod, dmod, 1);
            _perm_inv(P, P, n);
            v1mod = n_mulmod2_preinv(q, nmod_mat_entry(cmod, 0, P[n - 1]), p, Bmod->mod.ninv);
            v2mod = n_mulmod2_preinv(q, nmod_mat_entry(dmod, 0, P[n - 1]), p, Bmod->mod.ninv);*/

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
        flint_free(P);
        nmod_mat_clear(Btmod);
        /*nmod_mat_clear(Bmod);
        nmod_mat_clear(cmod);
        nmod_mat_clear(dmod);*/
    }
    else
    {
        fmpz_mat_det(d1, Bt);
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(Bt, j, n - 1), fmpz_mat_entry(d, 0, j));
        fmpz_mat_det(d2, Bt);
    }

    fmpz_clear(M);
    fmpz_mat_clear(dt);
    fmpz_mat_clear(Bt);
    fmpz_mat_clear(x);
    fmpq_mat_clear(x_q);
}

void fmpz_mat_hnf_pernet_stein(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong i, j;
    fmpz_t d1, d2, g, s, t;
    fmpz_mat_t c, d, B, C, H1, H2, H3, H4;

    fmpz_init(d1);
    fmpz_init(d2);
    fmpz_init(g);
    fmpz_init(s);
    fmpz_init(t);
    fmpz_mat_init(c, 1, A->c - 1);
    fmpz_mat_init(d, 1, A->c - 1);
    fmpz_mat_init(B, A->r - 2, A->c - 1);
    fmpz_mat_init(C, A->r - 1, A->c - 1);
    fmpz_mat_init(H1, A->r - 1, A->c - 1);
    fmpz_mat_init(H2, A->r - 1, A->c);
    fmpz_mat_init(H3, A->r, A->c);
    fmpz_mat_init(H4, A->r + 1, A->c);

    for (i = 0; i < A->r - 2; i++)
    {
        for (j = 0; j < A->c - 1; j++)
        {
            fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j));
            fmpz_set(fmpz_mat_entry(C, i, j), fmpz_mat_entry(A, i, j));
        }
        fmpz_set(fmpz_mat_entry(C, i, j), fmpz_mat_entry(A, i, j));
    }
    /*flint_printf("A:\n");
    fmpz_mat_print_pretty(A);
    flint_printf("\nB:\n");
    fmpz_mat_print_pretty(B);*/
    for (j = 0; j < A->c - 1; j++)
    {
        fmpz_set(fmpz_mat_entry(c, 0, j), fmpz_mat_entry(A, A->r - 2, j));
        fmpz_set(fmpz_mat_entry(d, 0, j), fmpz_mat_entry(A, A->r - 1, j));
    }
    double_det(d1, d2, B, c, d);
    fmpz_xgcd(g, s, t, d1, d2);
    for (j = 0; j < A->c - 1; j++)
    {
        fmpz_mul(fmpz_mat_entry(C, A->r - 2, j), s,
                fmpz_mat_entry(A, A->r - 2, j));
        fmpz_addmul(fmpz_mat_entry(C, A->r - 2, j), t,
                fmpz_mat_entry(A, A->r - 1, j));
    }
    if (!fmpz_is_zero(g)) /* chosen matrix invertible */
    {
        /* if g is too big, recurse */
        fmpz_abs(g, g);
        if (COEFF_IS_MPZ(*g) && 0)
        {
            fmpz_print(g);
            flint_printf("\n\n");
        }
        if (COEFF_IS_MPZ(*g) && C->r > 3)
        {
            fmpz_mat_hnf_pernet_stein(H1, C);
        }
        else /* use modulo determinant algorithm to compute HNF of C */
        {
            fmpz_mat_hnf_modular(H1, C, g);
        }
        fmpz_mat_clear(B);
        fmpz_mat_init(B, A->r - 1, A->c);
        for (j = 0; j < A->c; j++)
        {
            for (i = 0; i < A->r - 2; i++)
                fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j));
            fmpz_mul(fmpz_mat_entry(B, A->r - 2, j), s,
                    fmpz_mat_entry(A, A->r - 2, j));
            fmpz_addmul(fmpz_mat_entry(B, A->r - 2, j), t,
                    fmpz_mat_entry(A, A->r - 1, j));
        }
        add_column(H2, B, H1);
        add_row(H3, H2, A->rows[A->r - 2]);
        add_row(H4, H3, A->rows[A->r - 1]);
        /*flint_printf("\nd1/d2:\n");
        fmpz_print(d1);flint_printf("\n");
        fmpz_print(d2);
        flint_printf("\nB:\n");
        fmpz_mat_print_pretty(B);
        flint_printf("\nc:\n");
        fmpz_mat_print_pretty(c);
        flint_printf("\nd:\n");
        fmpz_mat_print_pretty(d);
        flint_printf("\nC:\n");
        fmpz_mat_print_pretty(C);
        flint_printf("\nH1:\n");
        fmpz_mat_print_pretty(H1);
        flint_printf("H2:\n");
        fmpz_mat_print_pretty(H2);
        flint_printf("H3:\n");
        fmpz_mat_print_pretty(H3);
        flint_printf("H4:\n");
        fmpz_mat_print_pretty(H4);*/
        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H4, i, j));
    }
    else
    {
        fmpz_mat_hnf_xgcd(H, A);
    }
    fmpz_clear(d2);
    fmpz_clear(d1);
    fmpz_clear(t);
    fmpz_clear(s);
    fmpz_clear(g);
    fmpz_mat_clear(H1);
    fmpz_mat_clear(H2);
    fmpz_mat_clear(H3);
    fmpz_mat_clear(H4);
    fmpz_mat_clear(C);
    fmpz_mat_clear(B);
    fmpz_mat_clear(c);
    fmpz_mat_clear(d);
}
