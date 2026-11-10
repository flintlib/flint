/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_vec.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* The evaluation of linear recurrences with sums and a scalar-geometric
   forcing against the exact rational iteration: random small (sparse)
   integer leaves of up to several thousand bits (several leaves per
   chunk, or one, and the transformed products), random precisions, one
   to three columns of the state, with and without a bottom function;
   the results must enclose the exact values, to about half the
   precision at least. */

typedef struct
{
    slong d, s, g;
    ulong seed;
    slong bits;
}
_rec_ctx;

/* a deterministic pseudo-random leaf from (seed, l) */
static void
_leaf(fmpz_mat_t T, fmpz_t phi, fmpz_t q, slong l, const void * vctx)
{
    const _rec_ctx * c = (const _rec_ctx *) vctx;
    flint_rand_t st;
    slong i, j;

    flint_rand_init(st);
    flint_rand_set_seed(st, c->seed, (ulong) l + 1);
    for (i = 0; i < c->d + c->s; i++)
        for (j = 0; j < c->d + c->g; j++)
            if (n_randint(st, 3) != 0)
                fmpz_randtest(fmpz_mat_entry(T, i, j), st, c->bits);
    if (c->g > 0)
        fmpz_randtest(phi, st, c->bits);
    /* |q| larger than the entries, so the products stay bounded */
    fmpz_randtest_not_zero(q, st, c->bits);
    fmpz_mul_2exp(q, q, 4 + FLINT_BIT_COUNT(c->d + c->g));
    fmpz_add_ui(q, q, 1);
    flint_rand_clear(st);
}

/* the same nodes by leaves, as a bottom function */
static void
_bottom(fmpz_mat_t T, fmpz_t phi, fmpz_t q, slong a, slong b, const void * vctx)
{
    const _rec_ctx * c = (const _rec_ctx *) vctx;
    fmpq_mat_t P, L;
    fmpq_t ph;
    fmpz_mat_t Tl;
    fmpz_t pl, ql;
    slong i, j, l, r = c->d + c->g + c->s;

    /* the node as a rational (d + g + s) x (d + g) block matrix acting on
       (X, f v) with the sums appended */
    fmpq_mat_init(P, r, c->d + c->g);
    fmpq_mat_init(L, r, r);
    fmpq_init(ph);
    fmpz_mat_init(Tl, c->d + c->s, c->d + c->g);
    fmpz_init(pl);
    fmpz_init(ql);
    for (i = 0; i < c->d + c->g; i++)
        fmpq_one(fmpq_mat_entry(P, i, i));
    for (l = a; l < b; l++)
    {
        fmpq_mat_t Q;
        fmpz_mat_zero(Tl);
        fmpz_zero(pl);
        _leaf(Tl, pl, ql, l, c);
        fmpq_mat_zero(L);
        for (i = 0; i < c->d; i++)
            for (j = 0; j < c->d + c->g; j++)
                fmpq_set_fmpz_frac(fmpq_mat_entry(L, i, j), fmpz_mat_entry(Tl, i, j), ql);
        for (i = 0; i < c->g; i++)
            fmpq_set_fmpz_frac(fmpq_mat_entry(L, c->d + i, c->d + i), pl, ql);
        for (i = 0; i < c->s; i++)
        {
            for (j = 0; j < c->d + c->g; j++)
                fmpq_set_fmpz_frac(fmpq_mat_entry(L, c->d + c->g + i, j), fmpz_mat_entry(Tl, c->d + i, j), ql);
            fmpq_one(fmpq_mat_entry(L, c->d + c->g + i, c->d + c->g + i));
        }
        /* P = L [P with the sums as extra columns ...]: keep P as
           (d + g + s) x (d + g + s) implicitly by applying to the
           augmented identity */
        fmpq_mat_init(Q, r, c->d + c->g);
        {
            fmpq_mat_t Pa;
            fmpq_mat_init(Pa, r, c->d + c->g);
            fmpq_mat_set(Pa, P);
            fmpq_mat_mul(Q, L, Pa);
            fmpq_mat_clear(Pa);
        }
        fmpq_mat_swap(P, Q);
        fmpq_mat_clear(Q);
    }
    /* back to the integer format over a common denominator */
    fmpz_one(q);
    for (i = 0; i < r; i++)
        for (j = 0; j < c->d + c->g; j++)
            fmpz_lcm(q, q, fmpq_denref(fmpq_mat_entry(P, i, j)));
    for (i = 0; i < c->d; i++)
        for (j = 0; j < c->d + c->g; j++)
        {
            fmpz_divexact(ql, q, fmpq_denref(fmpq_mat_entry(P, i, j)));
            fmpz_mul(fmpz_mat_entry(T, i, j), fmpq_numref(fmpq_mat_entry(P, i, j)), ql);
        }
    for (i = 0; i < c->s; i++)
        for (j = 0; j < c->d + c->g; j++)
        {
            fmpz_divexact(ql, q, fmpq_denref(fmpq_mat_entry(P, c->d + c->g + i, j)));
            fmpz_mul(fmpz_mat_entry(T, c->d + i, j), fmpq_numref(fmpq_mat_entry(P, c->d + c->g + i, j)), ql);
        }
    if (c->g > 0)
    {
        /* f_b / f_a */
        fmpz_divexact(ql, q, fmpq_denref(fmpq_mat_entry(P, c->d, c->d)));
        fmpz_mul(phi, fmpq_numref(fmpq_mat_entry(P, c->d, c->d)), ql);
    }
    fmpq_mat_clear(P);
    fmpq_mat_clear(L);
    fmpq_clear(ph);
    fmpz_mat_clear(Tl);
    fmpz_clear(pl);
    fmpz_clear(ql);
}

TEST_FUNCTION_START(mp_real_linrec, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        _rec_ctx c;
        _mp_real_linrec_struct rec;
        slong a, b, l, i, j, k, W, n, r;
        fmpq * x, * y, * v, * S;
        fmpq_t f, t, acc;
        mp_real_struct * xr, * vr, * sr;
        mp_real_t fr;
        fmpz_mat_t T;
        fmpz_t phi, q;

        c.d = 1 + n_randint(state, 4);
        c.s = n_randint(state, 4);
        c.g = n_randint(state, 3);
        c.seed = n_randlimb(state);
        c.bits = 2 + n_randint(state, (iter % 8 == 0) ? 3000 : 100);
        r = (c.g == 0) ? 1 + n_randint(state, 3) : 1;
        rec.d = c.d;
        rec.s = c.s;
        rec.g = c.g;
        rec.leaf = _leaf;
        rec.ctx = &c;
        if (n_randint(state, 2))
        {
            rec.bottom = _bottom;
            rec.bottom_len = 1 + n_randint(state, 5);
        }
        else
        {
            rec.bottom = NULL;
            rec.bottom_len = 0;
        }

        a = n_randint(state, 10);
        b = a + 1 + n_randint(state, (iter % 8 == 0) ? 60 : 20);
        W = 64 + n_randint(state, 500);
        /* now and then chunks beyond the cutoff for the transformed
           products (12000 bits) */
        if (iter % 16 == 7)
        {
            c.bits = 4000 + n_randint(state, 4000);
            b = a + 4 + n_randint(state, 12);
            W = 60000 + n_randint(state, 40000);
        }
        n = mp_real_prec_bits(W) + 4;

        x = _fmpq_vec_init(c.d * r);
        y = _fmpq_vec_init(c.d);
        v = _fmpq_vec_init(c.g);
        S = _fmpq_vec_init(c.s * r);
        fmpq_init(f);
        fmpq_init(t);
        fmpq_init(acc);
        xr = flint_malloc(sizeof(mp_real_struct) * c.d * r);
        vr = flint_malloc(sizeof(mp_real_struct) * FLINT_MAX(c.g, 1));
        sr = flint_malloc(sizeof(mp_real_struct) * FLINT_MAX(c.s * r, 1));
        for (i = 0; i < c.d * r; i++)
            mp_real_init(xr + i);
        for (i = 0; i < c.g; i++)
            mp_real_init(vr + i);
        for (i = 0; i < c.s * r; i++)
            mp_real_init(sr + i);
        mp_real_init(fr);

        /* dyadic inputs (including the initial sums), exact as balls */
        for (i = 0; i < c.d * r + c.s * r; i++)
        {
            fmpq * e = (i < c.d * r) ? x + i : S + (i - c.d * r);
            mp_real_struct * z = (i < c.d * r) ? xr + i : sr + (i - c.d * r);

            fmpz_randtest(fmpq_numref(e), state, 30);
            fmpz_one(fmpq_denref(e));
            fmpz_mul_2exp(fmpq_denref(e), fmpq_denref(e), n_randint(state, 20));
            fmpq_canonicalise(e);
            mp_real_set_fmpz(z, fmpq_numref(e));
            mp_real_mul_2exp_si(z, z, -(slong) fmpz_val2(fmpq_denref(e)));
        }
        for (i = 0; i < c.g; i++)
        {
            fmpz_randtest(fmpq_numref(v + i), state, 30);
            mp_real_set_fmpz(vr + i, fmpq_numref(v + i));
        }
        fmpz_randtest(fmpq_numref(f), state, 30);
        mp_real_set_fmpz(fr, fmpq_numref(f));

        /* the exact iteration, column by column */
        fmpz_mat_init(T, c.d + c.s, c.d + c.g);
        fmpz_init(phi);
        fmpz_init(q);
        for (l = a; l < b; l++)
        {
            fmpz_mat_zero(T);
            fmpz_zero(phi);
            _leaf(T, phi, q, l, &c);
            for (k = 0; k < r; k++)
            {
                for (i = 0; i < c.d + c.s; i++)
                {
                    fmpq_zero(acc);
                    for (j = 0; j < c.d; j++)
                    {
                        fmpq_mul_fmpz(t, x + j * r + k, fmpz_mat_entry(T, i, j));
                        fmpq_add(acc, acc, t);
                    }
                    for (j = 0; j < c.g; j++)
                    {
                        fmpq_mul_fmpz(t, v + j, fmpz_mat_entry(T, i, c.d + j));
                        fmpq_mul(t, t, f);
                        fmpq_add(acc, acc, t);
                    }
                    fmpq_div_fmpz(acc, acc, q);
                    if (i < c.d)
                        fmpq_set(y + i, acc);
                    else
                        fmpq_add(S + (i - c.d) * r + k, S + (i - c.d) * r + k, acc);
                }
                for (i = 0; i < c.d; i++)
                    fmpq_swap(x + i * r + k, y + i);
            }
            if (c.g > 0)
            {
                fmpq_mul_fmpz(f, f, phi);
                fmpq_div_fmpz(f, f, q);
            }
        }

        _mp_real_linrec_eval(xr, sr, (c.g > 0) ? fr : NULL, &rec, a, b,
            (c.g > 0) ? vr : NULL, r, n);

        /* enclosures, and no more than about half the bits lost */
        for (i = 0; i < c.d * r + c.s * r + (c.g > 0); i++)
        {
            arb_t A, E;
            mag_t m;
            const fmpq * ex = (i < c.d * r) ? x + i : (i < c.d * r + c.s * r) ? S + (i - c.d * r) : f;
            const mp_real_struct * z = (i < c.d * r) ? xr + i : (i < c.d * r + c.s * r) ? sr + (i - c.d * r) : fr;
            arb_init(A);
            arb_init(E);
            mag_init(m);
            mp_real_get_arb(A, z);
            arb_set_fmpq(E, ex, W + 200);
            arb_get_mag(m, E);
            mag_one(arb_radref(E));
            mag_max(m, m, arb_radref(E));
            mag_mul_2exp_si(m, m, -W / 2);
            if (!arb_contains_fmpq(A, ex) || mag_cmp(arb_radref(A), m) > 0)
            {
                flint_printf("FAIL: d = %wd, s = %wd, g = %wd, r = %wd, bits = %wd, [%wd, %wd), W = %wd, bottom = %d, entry %wd\n",
                    c.d, c.s, c.g, r, c.bits, a, b, W, rec.bottom != NULL, i);
                flint_printf("res = "); arb_printn(A, 30, 0); flint_printf("\n");
                arb_set_fmpq(E, ex, W + 200);
                flint_printf("exact = "); arb_printn(E, 30, 0); flint_printf("\n");
                flint_abort();
            }
            arb_clear(A);
            arb_clear(E);
            mag_clear(m);
        }

        fmpz_mat_clear(T);
        fmpz_clear(phi);
        fmpz_clear(q);
        _fmpq_vec_clear(x, c.d * r);
        _fmpq_vec_clear(y, c.d);
        _fmpq_vec_clear(v, c.g);
        _fmpq_vec_clear(S, c.s * r);
        fmpq_clear(f);
        fmpq_clear(t);
        fmpq_clear(acc);
        for (i = 0; i < c.d * r; i++)
            mp_real_clear(xr + i);
        for (i = 0; i < c.g; i++)
            mp_real_clear(vr + i);
        for (i = 0; i < c.s * r; i++)
            mp_real_clear(sr + i);
        mp_real_clear(fr);
        flint_free(xr);
        flint_free(vr);
        flint_free(sr);
    }

    TEST_FUNCTION_END(state);
}
