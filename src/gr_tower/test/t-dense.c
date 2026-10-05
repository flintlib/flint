/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "gr_special.h"
#include "fmpz_poly.h"
#include "qqbar.h"
#include "acb.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/impl.h"

#define DENSE_OK(expr) do { int __st = (expr); if (__st != GR_SUCCESS) { flint_printf("FAIL: %s returned %d (line %d)\n", #expr, __st, __LINE__); flint_abort(); } } while (0)

/*
    Dense arithmetic of lazy fields over number fields with univariate
    moduli (Q(zeta_n), quadratic fields, products of such): products,
    inverses, quotients, polynomial and matrix products agree with the
    generic algorithms, with entries from different towers and rational
    entries mixed in.
*/

/* a random element: a small combination of powers of the generators
   (roots of unity and square roots, chosen by kind) */
static void
_rand_elem(gr_ptr x, gr_srcptr const * gens, slong ngens, flint_rand_t state, gr_ctx_t ctx)
{
    gr_ptr t = gr_heap_init(ctx);
    slong i, terms = n_randint(state, 5);

    DENSE_OK(gr_set_si(x, (slong) n_randint(state, 7) - 3, ctx));
    if (n_randint(state, 4) == 0)
        DENSE_OK(gr_div_ui(x, x, 1 + n_randint(state, 5), ctx));

    for (i = 0; i < terms; i++)
    {
        DENSE_OK(gr_pow_ui(t, gens[n_randint(state, ngens)], n_randint(state, 12), ctx));
        if (n_randint(state, 2))
            DENSE_OK(gr_mul(t, t, gens[n_randint(state, ngens)], ctx));
        DENSE_OK(gr_mul_si(t, t, (slong) n_randint(state, 11) - 5, ctx));
        if (n_randint(state, 5) == 0)
            DENSE_OK(gr_div_ui(t, t, 1 + n_randint(state, 6), ctx));
        DENSE_OK(gr_add(x, x, t, ctx));
    }

    gr_heap_clear(t, ctx);
}

TEST_FUNCTION_START(gr_tower_dense, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K;
        gr_ptr gens[3];
        slong ngens, i, j;
        ulong orders[] = {3, 4, 5, 7, 8, 9, 12, 15, 16};

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        if (n_randint(state, 2))
            gr_tower_lazy_ctx_set_gen_flags(K, GR_TOWER_GENS_SPLIT_IMAGINARY);
        /* (the inverse algorithms: automatic, modular only, without the
           modular one) */
        GR_MUST_SUCCEED(gr_tower_lazy_ctx_set_option(K, GR_TOWER_OPT_INV_DENSE_ALG, n_randint(state, 3)));

        /* generators: roots of unity exp(2 pi i / n) and square roots */
        ngens = 1 + n_randint(state, 3);
        for (i = 0; i < ngens; i++)
        {
            gens[i] = gr_heap_init(K);
            if (n_randint(state, 3) == 0)
            {
                DENSE_OK(gr_set_si(gens[i], (n_randint(state, 2) ? 1 : -1) * (slong) (2 + n_randint(state, 20)), K));
                DENSE_OK(gr_sqrt(gens[i], gens[i], K));
            }
            else
            {
                gr_ptr t = gr_heap_init(K);
                DENSE_OK(gr_pi(gens[i], K));
                DENSE_OK(gr_i(t, K));
                DENSE_OK(gr_mul(gens[i], gens[i], t, K));
                DENSE_OK(gr_mul_ui(gens[i], gens[i], 2, K));
                DENSE_OK(gr_div_ui(gens[i], gens[i], orders[n_randint(state, 9)], K));
                DENSE_OK(gr_exp(gens[i], gens[i], K));
                gr_heap_clear(t, K);
            }
        }

        /* elements: products, inverses and quotients */
        {
            gr_ptr x, y, z, w;
            GR_TMP_INIT4(x, y, z, w, K);
            for (j = 0; j < 5; j++)
            {
                _rand_elem(x, (gr_srcptr const *) gens, ngens, state, K);
                /* (the same divisor again, possibly after extending the
                   tower: the cached inverse) */
                if (j == 0 || n_randint(state, 2))
                    _rand_elem(y, (gr_srcptr const *) gens, ngens, state, K);

                DENSE_OK(gr_mul(z, x, y, K));
                /* the product against the numerical values */
                {
                    acb_t vx, vy, vz;
                    acb_init(vx); acb_init(vy); acb_init(vz);
                    DENSE_OK(gr_tower_lazy_get_acb(vx, x, 64, K));
                    DENSE_OK(gr_tower_lazy_get_acb(vy, y, 64, K));
                    DENSE_OK(gr_tower_lazy_get_acb(vz, z, 64, K));
                    acb_mul(vx, vx, vy, 64);
                    if (!acb_overlaps(vx, vz))
                    {
                        flint_printf("FAIL: product (value)\n");
                        gr_println(x, K); gr_println(y, K); gr_println(z, K);
                        flint_abort();
                    }
                    acb_clear(vx); acb_clear(vy); acb_clear(vz);
                }
                if (gr_is_zero(y, K) == T_FALSE)
                {
                    DENSE_OK(gr_div(w, z, y, K));
                    if (gr_equal(w, x, K) != T_TRUE)
                    {
                        flint_printf("FAIL: (x y) / y\n");
                        gr_println(x, K); gr_println(y, K); gr_println(w, K);
                        flint_abort();
                    }
                    DENSE_OK(gr_inv(w, y, K));
                    DENSE_OK(gr_mul(w, w, y, K));
                    if (gr_is_one(w, K) != T_TRUE)
                    {
                        flint_printf("FAIL: y / y\n");
                        gr_println(y, K); gr_println(w, K);
                        flint_abort();
                    }
                }
            }
            GR_TMP_CLEAR4(x, y, z, w, K);
        }

        /* polynomial products */
        {
            gr_poly_t f, g, h1, h2;
            slong len1 = 1 + n_randint(state, 8), len2 = 1 + n_randint(state, 8), n;
            gr_poly_init(f, K);
            gr_poly_init(g, K);
            gr_poly_init(h1, K);
            gr_poly_init(h2, K);
            gr_poly_fit_length(f, len1, K);
            gr_poly_fit_length(g, len2, K);
            for (j = 0; j < len1; j++)
                _rand_elem(gr_poly_coeff_ptr(f, j, K), (gr_srcptr const *) gens, ngens, state, K);
            for (j = 0; j < len2; j++)
                _rand_elem(gr_poly_coeff_ptr(g, j, K), (gr_srcptr const *) gens, ngens, state, K);
            _gr_poly_set_length(f, len1, K);
            _gr_poly_set_length(g, len2, K);
            _gr_poly_normalise(f, K);
            _gr_poly_normalise(g, K);

            n = 1 + n_randint(state, len1 + len2);
            DENSE_OK(gr_poly_mullow(h1, f, g, n, K));
            if (f->length > 0 && g->length > 0)
            {
                slong m = FLINT_MIN(n, f->length + g->length - 1);
                gr_poly_fit_length(h2, m, K);
                DENSE_OK(_gr_poly_mullow_generic(h2->coeffs, f->coeffs, f->length, g->coeffs, g->length, m, K));
                _gr_poly_set_length(h2, m, K);
                _gr_poly_normalise(h2, K);
            }
            if (gr_poly_equal(h1, h2, K) != T_TRUE)
            {
                flint_printf("FAIL: polynomial product\n");
                gr_poly_print(f, K); flint_printf("\n");
                gr_poly_print(g, K); flint_printf("\n");
                gr_poly_print(h1, K); flint_printf("\n");
                gr_poly_print(h2, K); flint_printf("\n");
                flint_abort();
            }

            /* aliasing */
            DENSE_OK(gr_poly_mullow(f, f, g, n, K));
            if (gr_poly_equal(f, h2, K) != T_TRUE)
            {
                flint_printf("FAIL: polynomial product (aliasing)\n");
                flint_abort();
            }

            gr_poly_clear(f, K);
            gr_poly_clear(g, K);
            gr_poly_clear(h1, K);
            gr_poly_clear(h2, K);
        }

        /* matrix products */
        {
            gr_mat_t A, B, C1, C2;
            slong r = 1 + n_randint(state, 5), s = 1 + n_randint(state, 5), c = 1 + n_randint(state, 5);
            slong k;
            gr_mat_init(A, r, s, K);
            gr_mat_init(B, s, c, K);
            gr_mat_init(C1, r, c, K);
            gr_mat_init(C2, r, c, K);
            for (j = 0; j < r; j++)
                for (k = 0; k < s; k++)
                    _rand_elem(gr_mat_entry_ptr(A, j, k, K), (gr_srcptr const *) gens, ngens, state, K);
            for (j = 0; j < s; j++)
                for (k = 0; k < c; k++)
                    _rand_elem(gr_mat_entry_ptr(B, j, k, K), (gr_srcptr const *) gens, ngens, state, K);
            DENSE_OK(gr_mat_mul(C1, A, B, K));
            DENSE_OK(gr_mat_mul_classical(C2, A, B, K));
            if (gr_mat_equal(C1, C2, K) != T_TRUE)
            {
                flint_printf("FAIL: matrix product\n");
                gr_mat_print(A, K); gr_mat_print(B, K);
                gr_mat_print(C1, K); gr_mat_print(C2, K);
                flint_abort();
            }
            if (s == c)
            {
                /* aliasing */
                DENSE_OK(gr_mat_mul(A, A, B, K));
                if (gr_mat_equal(A, C2, K) != T_TRUE)
                {
                    flint_printf("FAIL: matrix product (aliasing)\n");
                    flint_abort();
                }
            }
            gr_mat_clear(A, K);
            gr_mat_clear(B, K);
            gr_mat_clear(C1, K);
            gr_mat_clear(C2, K);
        }

        /* determinants and linear systems (algorithm selection) */
        {
            gr_mat_t A, B, X, AX;
            gr_ptr d1, d2;
            slong n = 1 + n_randint(state, 7), k;
            gr_mat_init(A, n, n, K);
            gr_mat_init(B, n, 1 + n_randint(state, 2), K);
            gr_mat_init(X, n, gr_mat_ncols(B, K), K);
            gr_mat_init(AX, n, gr_mat_ncols(B, K), K);
            GR_TMP_INIT2(d1, d2, K);
            for (j = 0; j < n; j++)
            {
                for (k = 0; k < n; k++)
                    _rand_elem(gr_mat_entry_ptr(A, j, k, K), (gr_srcptr const *) gens, ngens, state, K);
                for (k = 0; k < gr_mat_ncols(B, K); k++)
                    _rand_elem(gr_mat_entry_ptr(B, j, k, K), (gr_srcptr const *) gens, ngens, state, K);
            }
            DENSE_OK(gr_mat_det(d1, A, K));
            DENSE_OK(gr_mat_det_berkowitz(d2, A, K));
            if (gr_equal(d1, d2, K) != T_TRUE)
            {
                flint_printf("FAIL: determinant\n");
                gr_mat_print(A, K); gr_println(d1, K); gr_println(d2, K);
                flint_abort();
            }
            if (gr_is_zero(d1, K) == T_FALSE)
            {
                DENSE_OK(gr_mat_nonsingular_solve(X, A, B, K));
                DENSE_OK(gr_mat_mul_classical(AX, A, X, K));
                if (gr_mat_equal(AX, B, K) != T_TRUE)
                {
                    flint_printf("FAIL: linear system\n");
                    gr_mat_print(A, K); gr_mat_print(B, K); gr_mat_print(X, K);
                    flint_abort();
                }
            }
            GR_TMP_CLEAR2(d1, d2, K);
            gr_mat_clear(A, K);
            gr_mat_clear(B, K);
            gr_mat_clear(X, K);
            gr_mat_clear(AX, K);
        }

        /* division by the same element before and after extending the
           tower (the cached inverse must not survive the extension) */
        {
            gr_ptr x, y, z, w;
            GR_TMP_INIT4(x, y, z, w, K);
            _rand_elem(y, (gr_srcptr const *) gens, ngens, state, K);
            if (gr_is_zero(y, K) == T_FALSE)
            {
                for (j = 0; j < 2; j++)
                {
                    _rand_elem(x, (gr_srcptr const *) gens, ngens, state, K);
                    if (j == 1)
                    {
                        DENSE_OK(gr_set_ui(z, 2 + n_randint(state, 30), K));
                        DENSE_OK(gr_sqrt(z, z, K));
                        DENSE_OK(gr_add(x, x, z, K));
                    }
                    DENSE_OK(gr_mul(z, x, y, K));
                    DENSE_OK(gr_div(w, z, y, K));
                    if (gr_equal(w, x, K) != T_TRUE)
                    {
                        flint_printf("FAIL: (x y) / y after extension\n");
                        gr_println(x, K); gr_println(y, K); gr_println(w, K);
                        flint_abort();
                    }
                }
            }
            GR_TMP_CLEAR4(x, y, z, w, K);
        }

        for (i = 0; i < ngens; i++)
            gr_heap_clear(gens[i], K);
        gr_ctx_clear(K);
    }

    /* inverses in a large tower: zeta_16, zeta_9, zeta_5 and a square
       root (degree 384) */
    for (iter = 0; iter < flint_test_multiplier(); iter++)
    {
        gr_ctx_t K;
        gr_ptr g, x, y, z, w;
        ulong orders[3] = {16, 9, 5};
        slong j;

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(g, x, y, z, w, K);
        DENSE_OK(gr_set_ui(x, 1 + n_randint(state, 5), K));
        DENSE_OK(gr_set_ui(y, 1, K));
        for (j = 0; j < 4; j++)
        {
            if (j < 3)
            {
                DENSE_OK(gr_pi(g, K));
                DENSE_OK(gr_i(z, K));
                DENSE_OK(gr_mul(g, g, z, K));
                DENSE_OK(gr_mul_ui(g, g, 2, K));
                DENSE_OK(gr_div_ui(g, g, orders[j], K));
                DENSE_OK(gr_exp(g, g, K));
            }
            else
            {
                DENSE_OK(gr_set_ui(g, 7, K));
                DENSE_OK(gr_sqrt(g, g, K));
            }
            DENSE_OK(gr_add_si(z, g, (slong) n_randint(state, 5) - 2, K));
            DENSE_OK(gr_mul(x, x, z, K));
            DENSE_OK(gr_add_si(z, g, (slong) n_randint(state, 7) - 3, K));
            DENSE_OK(gr_add(y, y, z, K));
            DENSE_OK(gr_mul(y, y, g, K));
            DENSE_OK(gr_add_si(y, y, 1, K));
        }
        DENSE_OK(gr_add(x, x, y, K));
        if (gr_is_zero(y, K) == T_FALSE)
        {
            DENSE_OK(gr_mul(z, x, y, K));
            DENSE_OK(gr_div(w, z, y, K));
            if (gr_equal(w, x, K) != T_TRUE)
            {
                flint_printf("FAIL: large tower\n");
                gr_println(x, K); gr_println(y, K); gr_println(w, K);
                flint_abort();
            }
            DENSE_OK(gr_inv(w, y, K));
            DENSE_OK(gr_mul(w, w, y, K));
            if (gr_is_one(w, K) != T_TRUE)
            {
                flint_printf("FAIL: large tower (inverse)\n");
                flint_abort();
            }
        }
        GR_TMP_CLEAR5(g, x, y, z, w, K);
        gr_ctx_clear(K);
    }

    /* square roots of negative rational numbers: a single generator
       sqrt(-A) by default, i sqrt(A) with GR_TOWER_GENS_SPLIT_IMAGINARY;
       the two agree, and products like sqrt(-a) sqrt(-b) = -sqrt(a b)
       are recognized */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K;
        gr_ptr x, y, z, w;
        fmpq_t a, b;

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT4(x, y, z, w, K);
        fmpq_init(a);
        fmpq_init(b);
        fmpq_set_si(a, 1 + n_randint(state, 200), 1 + n_randint(state, 30));
        fmpq_set_si(b, 1 + n_randint(state, 200), 1 + n_randint(state, 30));

        /* x = sqrt(-a), y = sqrt(-b), with the flags in some order */
        if (n_randint(state, 2))
            gr_tower_lazy_ctx_set_gen_flags(K, GR_TOWER_GENS_SPLIT_IMAGINARY);
        DENSE_OK(gr_set_fmpq(x, a, K));
        DENSE_OK(gr_neg(x, x, K));
        DENSE_OK(gr_sqrt(x, x, K));
        gr_tower_lazy_ctx_set_gen_flags(K, n_randint(state, 2) ? GR_TOWER_GENS_SPLIT_IMAGINARY : 0);
        DENSE_OK(gr_set_fmpq(y, b, K));
        DENSE_OK(gr_neg(y, y, K));
        DENSE_OK(gr_sqrt(y, y, K));

        /* x^2 = -a */
        DENSE_OK(gr_sqr(z, x, K));
        DENSE_OK(gr_set_fmpq(w, a, K));
        DENSE_OK(gr_neg(w, w, K));
        if (gr_equal(z, w, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(-a)^2\n");
            flint_abort();
        }

        /* x y = -sqrt(a b) */
        DENSE_OK(gr_mul(z, x, y, K));
        fmpq_mul(b, a, b);
        DENSE_OK(gr_set_fmpq(w, b, K));
        DENSE_OK(gr_sqrt(w, w, K));
        DENSE_OK(gr_neg(w, w, K));
        if (gr_equal(z, w, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(-a) sqrt(-b)\n");
            gr_println(x, K); gr_println(y, K); gr_println(z, K); gr_println(w, K);
            flint_abort();
        }

        /* the imaginary part is sqrt(a); the real part is zero */
        DENSE_OK(gr_im(z, x, K));
        DENSE_OK(gr_set_fmpq(w, a, K));
        DENSE_OK(gr_sqrt(w, w, K));
        if (gr_equal(z, w, K) != T_TRUE)
        {
            flint_printf("FAIL: im(sqrt(-a))\n");
            flint_abort();
        }
        DENSE_OK(gr_re(z, x, K));
        if (gr_is_zero(z, K) != T_TRUE)
        {
            flint_printf("FAIL: re(sqrt(-a))\n");
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        GR_TMP_CLEAR4(x, y, z, w, K);
        gr_ctx_clear(K);
    }

    /* a large dense modulus (fast division): a root of the Eisenstein
       polynomial x^d + 2 (1 + 3 x + ... ) of degree 64 .. 80 */
    {
        gr_ctx_t K;
        fmpz_poly_t f;
        qqbar_ptr rts;
        gr_ptr a, x, y, z, w;
        slong d = 64 + n_randint(state, 17), j, k;

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        fmpz_poly_init(f);
        fmpz_poly_set_coeff_ui(f, d, 1);
        for (j = 0; j < d; j++)
            fmpz_poly_set_coeff_ui(f, j, 2 * (1 + n_randint(state, 5)));
        fmpz_poly_set_coeff_ui(f, 0, 2);
        rts = _qqbar_vec_init(d);
        qqbar_roots_fmpz_poly(rts, f, QQBAR_ROOTS_IRREDUCIBLE);

        GR_TMP_INIT5(a, x, y, z, w, K);
        {
            gr_ctx_t QQbar;
            gr_ctx_init_complex_qqbar(QQbar);
            DENSE_OK(gr_set_other(a, rts + 0, QQbar, K));
            gr_ctx_clear(QQbar);
        }

        for (k = 0; k < 3; k++)
        {
            DENSE_OK(gr_zero(x, K));
            DENSE_OK(gr_zero(y, K));
            for (j = d - 1; j >= 0; j--)
            {
                DENSE_OK(gr_mul(x, x, a, K));
                DENSE_OK(gr_add_si(x, x, (slong) n_randint(state, 100) - 50, K));
                if (j < 10)
                {
                    DENSE_OK(gr_mul(y, y, a, K));
                    DENSE_OK(gr_add_si(y, y, (slong) n_randint(state, 100) - 50, K));
                }
            }
            DENSE_OK(gr_mul(z, x, y, K));
            DENSE_OK(gr_div(w, z, x, K));
            if (gr_equal(w, y, K) != T_TRUE)
            {
                flint_printf("FAIL: dense modulus of degree %wd\n", d);
                flint_abort();
            }
            /* against the evaluation: z = x y numerically */
            {
                acb_t zx, zy, zz;
                acb_init(zx); acb_init(zy); acb_init(zz);
                DENSE_OK(gr_tower_lazy_get_acb(zx, x, 128, K));
                DENSE_OK(gr_tower_lazy_get_acb(zy, y, 128, K));
                DENSE_OK(gr_tower_lazy_get_acb(zz, z, 128, K));
                acb_mul(zx, zx, zy, 128);
                if (!acb_overlaps(zx, zz))
                {
                    flint_printf("FAIL: dense modulus of degree %wd (value)\n", d);
                    flint_abort();
                }
                acb_clear(zx); acb_clear(zy); acb_clear(zz);
            }
        }

        GR_TMP_CLEAR5(a, x, y, z, w, K);
        _qqbar_vec_clear(rts, d);
        fmpz_poly_clear(f);
        gr_ctx_clear(K);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
