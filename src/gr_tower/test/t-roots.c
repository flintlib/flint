/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

/* Structured roots in the lazy field: roots of unity and roots of
   integers as products of roots of prime power orders, and splitting
   towers for the roots of irreducible polynomials. */

/* exp(2 pi i p / q) */
static void
_root_of_unity(gr_ptr res, slong p, ulong q, gr_ctx_t K)
{
    gr_ptr t;
    t = gr_heap_init(K);
    GR_MUST_SUCCEED(gr_pi(res, K));
    GR_MUST_SUCCEED(gr_i(t, K));
    GR_MUST_SUCCEED(gr_mul(res, res, t, K));
    GR_MUST_SUCCEED(gr_mul_si(res, res, 2 * p, K));
    GR_MUST_SUCCEED(gr_div_ui(res, res, q, K));
    GR_MUST_SUCCEED(gr_exp(res, res, K));
    gr_heap_clear(t, K);
}

static void
_check_equal(gr_srcptr x, gr_srcptr y, const char * what, gr_ctx_t K)
{
    if (gr_equal(x, y, K) != T_TRUE)
    {
        flint_printf("FAIL: %s\n", what);
        gr_println(x, K);
        gr_println(y, K);
        flint_abort();
    }
}

TEST_FUNCTION_START(gr_tower_roots, state)
{
    gr_ctx_t QQ, K;
    gr_ptr x, y, z, t;
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
    x = gr_heap_init(K);
    y = gr_heap_init(K);
    z = gr_heap_init(K);
    t = gr_heap_init(K);

    /* roots of unity multiply by adding the fractions */
    for (iter = 0; iter < 100; iter++)
    {
        ulong q1 = 1 + n_randint(state, 30), q2 = 1 + n_randint(state, 30);
        slong p1 = n_randint(state, 60) - 30, p2 = n_randint(state, 60) - 30;
        fmpq_t a, b;

        fmpq_init(a);
        fmpq_init(b);
        fmpq_set_si(a, p1, q1);
        fmpq_set_si(b, p2, q2);
        fmpq_add(a, a, b);

        _root_of_unity(x, p1, q1, K);
        _root_of_unity(y, p2, q2, K);
        GR_MUST_SUCCEED(gr_mul(z, x, y, K));
        _root_of_unity(t, fmpz_get_si(fmpq_numref(a)), fmpz_get_ui(fmpq_denref(a)), K);
        _check_equal(z, t, "product of roots of unity", K);

        /* and the q-th power is one */
        GR_MUST_SUCCEED(gr_pow_ui(z, x, q1, K));
        GR_MUST_SUCCEED(gr_one(t, K));
        _check_equal(z, t, "root of unity to the power q", K);

        fmpq_clear(a);
        fmpq_clear(b);
    }

    /* the roots of unity of order 12 are powers of a common element:
       zeta_3 zeta_4 = zeta_12^7, and conjugation inverts them */
    _root_of_unity(x, 1, 3, K);
    _root_of_unity(y, 1, 4, K);
    GR_MUST_SUCCEED(gr_mul(z, x, y, K));
    _root_of_unity(t, 7, 12, K);
    _check_equal(z, t, "zeta_3 zeta_4", K);
    GR_MUST_SUCCEED(gr_conj(z, t, K));
    GR_MUST_SUCCEED(gr_inv(t, t, K));
    _check_equal(z, t, "conj of a root of unity", K);
    GR_MUST_SUCCEED(gr_abs(z, t, K));
    GR_MUST_SUCCEED(gr_one(t, K));
    _check_equal(z, t, "abs of a root of unity", K);

    /* i is exp(pi i / 2) */
    GR_MUST_SUCCEED(gr_i(x, K));
    _root_of_unity(y, 1, 4, K);
    _check_equal(x, y, "i", K);

    /* roots of integers: sqrt(12) = 2 sqrt(3), 6^(1/6) = 2^(1/6) 3^(1/6),
       2^(1/6) = sqrt(2) / 2^(1/3), sqrt(-2) = i sqrt(2) */
    GR_MUST_SUCCEED(gr_set_ui(x, 12, K));
    GR_MUST_SUCCEED(gr_sqrt(x, x, K));
    GR_MUST_SUCCEED(gr_set_ui(y, 3, K));
    GR_MUST_SUCCEED(gr_sqrt(y, y, K));
    GR_MUST_SUCCEED(gr_mul_ui(y, y, 2, K));
    _check_equal(x, y, "sqrt(12)", K);

    GR_MUST_SUCCEED(gr_set_ui(x, 6, K));
    GR_MUST_SUCCEED(gr_rsqrt(t, x, K));
    GR_MUST_SUCCEED(gr_sqrt(x, x, K));
    GR_MUST_SUCCEED(gr_mul(t, t, x, K));
    GR_MUST_SUCCEED(gr_one(z, K));
    _check_equal(t, z, "sqrt(6) rsqrt(6)", K);

    GR_MUST_SUCCEED(gr_set_ui(x, 6, K));
    GR_MUST_SUCCEED(gr_set_ui(t, 6, K));
    GR_MUST_SUCCEED(gr_tower_lazy_root_ui(x, x, 6, K));
    GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
    GR_MUST_SUCCEED(gr_tower_lazy_root_ui(y, y, 6, K));
    GR_MUST_SUCCEED(gr_set_ui(z, 3, K));
    GR_MUST_SUCCEED(gr_tower_lazy_root_ui(z, z, 6, K));
    GR_MUST_SUCCEED(gr_mul(y, y, z, K));
    _check_equal(x, y, "6^(1/6)", K);
    GR_MUST_SUCCEED(gr_pow_ui(y, x, 6, K));
    _check_equal(y, t, "(6^(1/6))^6", K);

    GR_MUST_SUCCEED(gr_set_ui(x, 2, K));
    GR_MUST_SUCCEED(gr_tower_lazy_root_ui(x, x, 6, K));
    GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
    GR_MUST_SUCCEED(gr_sqrt(y, y, K));
    GR_MUST_SUCCEED(gr_set_ui(z, 2, K));
    GR_MUST_SUCCEED(gr_tower_lazy_root_ui(z, z, 3, K));
    GR_MUST_SUCCEED(gr_div(y, y, z, K));
    _check_equal(x, y, "2^(1/6)", K);

    GR_MUST_SUCCEED(gr_set_si(x, -2, K));
    GR_MUST_SUCCEED(gr_sqrt(x, x, K));
    GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
    GR_MUST_SUCCEED(gr_sqrt(y, y, K));
    GR_MUST_SUCCEED(gr_i(z, K));
    GR_MUST_SUCCEED(gr_mul(y, y, z, K));
    _check_equal(x, y, "sqrt(-2)", K);

    /* roots of rational numbers: sqrt(9/8) = 3/4 sqrt(2) */
    {
        fmpq_t q;
        fmpq_init(q);
        fmpq_set_si(q, 9, 8);
        GR_MUST_SUCCEED(gr_set_fmpq(x, q, K));
        GR_MUST_SUCCEED(gr_sqrt(x, x, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(y, y, K));
        fmpq_set_si(q, 3, 4);
        GR_MUST_SUCCEED(gr_mul_fmpq(y, y, q, K));
        _check_equal(x, y, "sqrt(9/8)", K);
        fmpq_clear(q);
    }

    /* log of a Gaussian rational found through conjugation: log(1+i) =
       log(2)/2 + pi i / 4, even with i living inside Q(zeta_12) */
    _root_of_unity(x, 3, 12, K);
    GR_MUST_SUCCEED(gr_add_ui(x, x, 1, K));
    GR_MUST_SUCCEED(gr_log(x, x, K));
    GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
    GR_MUST_SUCCEED(gr_log(y, y, K));
    GR_MUST_SUCCEED(gr_div_ui(y, y, 2, K));
    GR_MUST_SUCCEED(gr_pi(z, K));
    GR_MUST_SUCCEED(gr_i(t, K));
    GR_MUST_SUCCEED(gr_mul(z, z, t, K));
    GR_MUST_SUCCEED(gr_div_ui(z, z, 4, K));
    GR_MUST_SUCCEED(gr_add(y, y, z, K));
    _check_equal(x, y, "log(1+i)", K);

    /* polynomial roots: x^n - 1 has the roots of unity, and the roots of
       an irreducible polynomial satisfy Vieta's formulas */
    {
        gr_poly_t f;
        gr_vec_t roots;
        fmpz_vec_t mult;
        slong i, n;

        gr_poly_init(f, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);

        for (n = 1; n <= 12; n++)
        {
            GR_MUST_SUCCEED(gr_poly_zero(f, K));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, n, 1, K));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, -1, K));
            GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));
            if (roots->length != n)
            {
                flint_printf("FAIL: number of roots of x^%wd - 1: %wd\n", n, roots->length);
                flint_abort();
            }
            for (i = 0; i < n; i++)
            {
                GR_MUST_SUCCEED(gr_pow_ui(x, gr_vec_entry_ptr(roots, i, K), n, K));
                GR_MUST_SUCCEED(gr_one(y, K));
                _check_equal(x, y, "root of x^n - 1", K);
            }
        }

        /* x^5 - x - 1, x^7 - 7 x + 3 and (x^3 - 2)(x^4 + x + 1)^2 */
        for (iter = 0; iter < 3; iter++)
        {
            gr_poly_t g;
            gr_poly_init(g, K);
            GR_MUST_SUCCEED(gr_poly_zero(f, K));
            if (iter == 0)
            {
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 5, 1, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 1, -1, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, -1, K));
            }
            else if (iter == 1)
            {
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 7, 1, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 1, -7, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, 3, K));
            }
            else
            {
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 3, 1, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, -2, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 4, 1, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 1, 1, K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 0, 1, K));
                GR_MUST_SUCCEED(gr_poly_mul(f, f, g, K));
                GR_MUST_SUCCEED(gr_poly_mul(f, f, g, K));
            }

            GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));

            /* the roots are roots, and prod (x - r_i)^m_i = f */
            GR_MUST_SUCCEED(gr_poly_one(g, K));
            for (i = 0; i < roots->length; i++)
            {
                gr_poly_t lin;
                slong j;
                GR_MUST_SUCCEED(gr_poly_evaluate(x, f, gr_vec_entry_ptr(roots, i, K), K));
                GR_MUST_SUCCEED(gr_zero(y, K));
                _check_equal(x, y, "evaluation at a root", K);
                gr_poly_init(lin, K);
                GR_MUST_SUCCEED(gr_poly_set_coeff_si(lin, 1, 1, K));
                GR_MUST_SUCCEED(gr_neg(x, gr_vec_entry_ptr(roots, i, K), K));
                GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(lin, 0, x, K));
                for (j = 0; j < fmpz_get_si(mult->entries + i); j++)
                    GR_MUST_SUCCEED(gr_poly_mul(g, g, lin, K));
                gr_poly_clear(lin, K);
            }
            if (gr_poly_equal(f, g, K) != T_TRUE)
            {
                flint_printf("FAIL: product of linear factors\n");
                gr_poly_print(f, K); flint_printf("\n");
                gr_poly_print(g, K); flint_printf("\n");
                flint_abort();
            }
            gr_poly_clear(g, K);
        }

        gr_poly_clear(f, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
    }

    /* a root of unity moved in front of a radical: sqrt(2) (1+i)/2 = exp(pi i/4),
       where x^2 - 2 becomes reducible over Q(zeta_8) (the proof of
       irreducibility must be dropped when the prefix grows) */
    {
        gr_ptr pi = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_pi(pi, K));
        GR_MUST_SUCCEED(gr_i(y, K));
        GR_MUST_SUCCEED(gr_add_ui(x, y, 1, K));
        GR_MUST_SUCCEED(gr_set_ui(z, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(z, z, K));
        GR_MUST_SUCCEED(gr_mul(x, z, x, K));
        GR_MUST_SUCCEED(gr_div_ui(x, x, 2, K));
        GR_MUST_SUCCEED(gr_mul(y, pi, y, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 4, K));
        GR_MUST_SUCCEED(gr_exp(y, y, K));
        _check_equal(x, y, "sqrt(2) (1+i)/2 = exp(pi i/4)", K);
        GR_MUST_SUCCEED(gr_mul(x, x, pi, K));
        GR_MUST_SUCCEED(gr_mul(y, y, pi, K));
        _check_equal(x, y, "sqrt(2) (1+i)/2 pi = exp(pi i/4) pi", K);

        GR_MUST_SUCCEED(gr_set_ui(x, 3, K));
        GR_MUST_SUCCEED(gr_sqrt(x, x, K));
        GR_MUST_SUCCEED(gr_i(y, K));
        GR_MUST_SUCCEED(gr_add(x, x, y, K));
        GR_MUST_SUCCEED(gr_div_ui(x, x, 2, K));
        GR_MUST_SUCCEED(gr_mul(y, pi, y, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 6, K));
        GR_MUST_SUCCEED(gr_exp(y, y, K));
        _check_equal(x, y, "(sqrt(3) + i)/2 = exp(pi i/6)", K);
        gr_heap_clear(pi, K);
    }

    /* exponentials of rational combinations of logarithms of rationals
       are products of radicals: exp(log(2)/3 + log(3)/2) = 2^(1/3) sqrt(3),
       exp((log(2i) - pi i/2)/2) = sqrt(2) (Calcium issue #24),
       exp(2 log 5 - log 3 + pi i/3) = 25/3 exp(pi i/3), exp(-log(8)/3) = 1/2 */
    {
        gr_ptr pi = gr_heap_init(K);
        fmpq_t q;
        fmpq_init(q);
        GR_MUST_SUCCEED(gr_pi(pi, K));

        GR_MUST_SUCCEED(gr_set_ui(x, 2, K));
        GR_MUST_SUCCEED(gr_log(x, x, K));
        GR_MUST_SUCCEED(gr_div_ui(x, x, 3, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 3, K));
        GR_MUST_SUCCEED(gr_log(y, y, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 2, K));
        GR_MUST_SUCCEED(gr_add(x, x, y, K));
        GR_MUST_SUCCEED(gr_exp(x, x, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
        fmpq_set_si(q, 1, 3);
        GR_MUST_SUCCEED(gr_pow_fmpq(y, y, q, K));
        GR_MUST_SUCCEED(gr_set_ui(z, 3, K));
        GR_MUST_SUCCEED(gr_sqrt(z, z, K));
        GR_MUST_SUCCEED(gr_mul(y, y, z, K));
        _check_equal(x, y, "exp(log(2)/3 + log(3)/2) = 2^(1/3) sqrt(3)", K);
        if (gr_tower_lazy_get_tower(NULL, x, K)->num_trans != 0)
        {
            flint_printf("FAIL: exp of logs not algebraic\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_i(y, K));
        GR_MUST_SUCCEED(gr_mul_ui(x, y, 2, K));
        GR_MUST_SUCCEED(gr_log(x, x, K));
        GR_MUST_SUCCEED(gr_mul(y, y, pi, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 2, K));
        GR_MUST_SUCCEED(gr_sub(x, x, y, K));
        GR_MUST_SUCCEED(gr_div_ui(x, x, 2, K));
        GR_MUST_SUCCEED(gr_exp(x, x, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(y, y, K));
        _check_equal(x, y, "exp((log(2i) - pi i/2)/2) = sqrt(2)", K);
        if (gr_tower_lazy_get_tower(NULL, x, K)->num_trans != 0)
        {
            flint_printf("FAIL: issue #24 value not algebraic\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_set_ui(x, 5, K));
        GR_MUST_SUCCEED(gr_log(x, x, K));
        GR_MUST_SUCCEED(gr_mul_ui(x, x, 2, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 3, K));
        GR_MUST_SUCCEED(gr_log(y, y, K));
        GR_MUST_SUCCEED(gr_sub(x, x, y, K));
        GR_MUST_SUCCEED(gr_i(y, K));
        GR_MUST_SUCCEED(gr_mul(y, y, pi, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 3, K));
        GR_MUST_SUCCEED(gr_add(x, x, y, K));
        GR_MUST_SUCCEED(gr_exp(x, x, K));
        GR_MUST_SUCCEED(gr_exp(y, y, K));
        GR_MUST_SUCCEED(gr_mul_ui(y, y, 25, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 3, K));
        _check_equal(x, y, "exp(2 log 5 - log 3 + pi i/3) = 25/3 exp(pi i/3)", K);

        GR_MUST_SUCCEED(gr_set_ui(x, 8, K));
        GR_MUST_SUCCEED(gr_log(x, x, K));
        GR_MUST_SUCCEED(gr_div_si(x, x, -3, K));
        GR_MUST_SUCCEED(gr_exp(x, x, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 1, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 2, K));
        _check_equal(x, y, "exp(-log(8)/3) = 1/2", K);

        fmpq_clear(q);
        gr_heap_clear(pi, K);
    }

    /* conjugation of a non-real root of an irreducible polynomial of
       degree 10 (Lehmer's polynomial), and re/im/abs of it */
    {
        gr_poly_t f;
        gr_vec_t roots;
        fmpz_vec_t mult;
        slong cs[11] = {1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1};
        slong i;

        gr_poly_init(f, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);
        for (i = 0; i <= 10; i++)
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, i, cs[i], K));
        GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));
        if (roots->length != 10)
        {
            flint_printf("FAIL: number of roots\n");
            flint_abort();
        }

        for (i = 0; i < 10; i++)
        {
            gr_srcptr r = gr_vec_entry_ptr(roots, i, K);
            GR_MUST_SUCCEED(gr_conj(x, r, K));
            GR_MUST_SUCCEED(gr_poly_evaluate(y, f, x, K));
            GR_MUST_SUCCEED(gr_zero(z, K));
            _check_equal(y, z, "f(conj(r)) = 0", K);
            GR_MUST_SUCCEED(gr_re(y, r, K));
            GR_MUST_SUCCEED(gr_im(z, r, K));
            GR_MUST_SUCCEED(gr_i(t, K));
            GR_MUST_SUCCEED(gr_mul(z, z, t, K));
            GR_MUST_SUCCEED(gr_add(y, y, z, K));
            _check_equal(y, r, "re(r) + i im(r) = r", K);
            GR_MUST_SUCCEED(gr_abs(y, r, K));
            GR_MUST_SUCCEED(gr_sqr(y, y, K));
            GR_MUST_SUCCEED(gr_mul(z, x, r, K));
            _check_equal(y, z, "|r|^2 = r conj(r)", K);
        }

        gr_poly_clear(f, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
    }

    /* roots linear in transcendental generators: (x - pi)(x - e)(x - pi e)
       (multivariate factorization over Q[pi, e]) */
    {
        gr_poly_t f, g;
        gr_vec_t roots;
        fmpz_vec_t mult;
        gr_ptr r[3];
        slong i, j;

        gr_poly_init(f, K);
        gr_poly_init(g, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);
        for (i = 0; i < 3; i++)
            r[i] = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_pi(r[0], K));
        GR_MUST_SUCCEED(gr_one(r[1], K));
        GR_MUST_SUCCEED(gr_exp(r[1], r[1], K));
        GR_MUST_SUCCEED(gr_mul(r[2], r[0], r[1], K));
        GR_MUST_SUCCEED(gr_poly_one(f, K));
        for (i = 0; i < 3; i++)
        {
            GR_MUST_SUCCEED(gr_poly_zero(g, K));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(g, 1, 1, K));
            GR_MUST_SUCCEED(gr_neg(x, r[i], K));
            GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(g, 0, x, K));
            GR_MUST_SUCCEED(gr_poly_mul(f, f, g, K));
        }
        GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));
        if (roots->length != 3)
        {
            flint_printf("FAIL: roots of (x - pi)(x - e)(x - pi e)\n");
            flint_abort();
        }
        for (i = 0; i < 3; i++)
        {
            slong hits = 0;
            for (j = 0; j < 3; j++)
                if (gr_equal(gr_vec_entry_srcptr(roots, j, K), r[i], K) == T_TRUE)
                    hits++;
            if (hits != 1)
            {
                flint_printf("FAIL: root %wd of (x - pi)(x - e)(x - pi e)\n", i);
                flint_abort();
            }
        }
        for (i = 0; i < 3; i++)
            gr_heap_clear(r[i], K);
        gr_poly_clear(f, K);
        gr_poly_clear(g, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
    }

    /* conjugates of roots of polynomials over transcendental generators
       (not absolute algebraic numbers): x^5 - pi x - 1 has three real
       roots; the conjugate of a nonreal root is the other nonreal root */
    {
        gr_poly_t f;
        gr_vec_t roots;
        fmpz_vec_t mult;
        slong i, j, nreal = 0;

        gr_poly_init(f, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);
        GR_MUST_SUCCEED(gr_pi(x, K));
        GR_MUST_SUCCEED(gr_neg(x, x, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 5, 1, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, 1, x, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, -1, K));
        GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));
        if (roots->length != 5)
        {
            flint_printf("FAIL: roots of x^5 - pi x - 1\n");
            flint_abort();
        }
        for (i = 0; i < 5; i++)
        {
            gr_srcptr r = gr_vec_entry_srcptr(roots, i, K);
            truth_t eq;
            GR_MUST_SUCCEED(gr_conj(y, r, K));
            eq = gr_equal(y, r, K);
            if (eq == T_TRUE)
                nreal++;
            else if (eq == T_FALSE)
            {
                slong hits = 0;
                for (j = 0; j < 5; j++)
                    if (j != i && gr_equal(y, gr_vec_entry_srcptr(roots, j, K), K) == T_TRUE)
                        hits++;
                if (hits != 1)
                {
                    flint_printf("FAIL: conj of a nonreal root of x^5 - pi x - 1\n");
                    flint_abort();
                }
            }
            else
            {
                flint_printf("FAIL: conj(r) = r undecided\n");
                flint_abort();
            }
            GR_MUST_SUCCEED(gr_conj(z, y, K));
            _check_equal(z, r, "conj(conj(r)) = r", K);
        }
        if (nreal != 3)
        {
            flint_printf("FAIL: three real roots of x^5 - pi x - 1 (%wd)\n", nreal);
            flint_abort();
        }
        gr_poly_clear(f, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
    }

    gr_heap_clear(x, K);
    gr_heap_clear(y, K);
    gr_heap_clear(z, K);
    gr_heap_clear(t, K);
    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
