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
#include "qqbar.h"
#include "gr_vec.h"
#include "acb.h"
#include "gr_poly.h"
#include "fmpz_vec.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

TEST_FUNCTION_START(gr_tower_lazy, state)
{
    gr_ctx_t QQ, K, QQbar;
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_complex_qqbar(QQbar);

    /* Generic ring tests, with and without express-on-merge */
    gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
    gr_test_ring(K, 5 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);
    gr_ctx_clear(K);

    gr_ctx_init_tower_lazy(K, QQ, 0);
    gr_test_ring(K, 5 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);
    gr_ctx_clear(K);

    /* Square roots: sqrt(2) sqrt(3) == sqrt(6), sqrt(4) == 2, sqrt(8) == 2 sqrt(2),
       and the towers involved. */
    {
        gr_ptr a, b, c, d, e;
        gr_tower_struct * T;
        slong level;

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(a, b, c, d, e, K);

        GR_MUST_SUCCEED(gr_set_ui(a, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(a, a, K));
        GR_MUST_SUCCEED(gr_set_ui(b, 3, K));
        GR_MUST_SUCCEED(gr_sqrt(b, b, K));
        GR_MUST_SUCCEED(gr_set_ui(c, 6, K));
        GR_MUST_SUCCEED(gr_sqrt(c, c, K));
        GR_MUST_SUCCEED(gr_mul(d, a, b, K));

        if (gr_equal(c, d, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(2) sqrt(3) != sqrt(6)\n");
            gr_println(c, K); gr_println(d, K);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_set_ui(e, 4, K));
        GR_MUST_SUCCEED(gr_sqrt(e, e, K));
        gr_tower_lazy_get_tower(&level, e, K);
        if (level != 0 || gr_sub_ui(e, e, 2, K) != GR_SUCCESS || gr_is_zero(e, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(4)\n");
            gr_println(e, K);
            flint_abort();
        }

        /* sqrt(8) computed from an element in the tower of sqrt(2) should be expressed */
        GR_MUST_SUCCEED(gr_mul_ui(e, a, 4, K));    /* 4 sqrt(2) */
        GR_MUST_SUCCEED(gr_sqr(e, e, K));           /* 32 */
        GR_MUST_SUCCEED(gr_sqrt(e, e, K));          /* sqrt(32) = 4 sqrt(2), a rational so a new tower */
        GR_MUST_SUCCEED(gr_mul_ui(d, a, 4, K));
        if (gr_equal(e, d, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(32) != 4 sqrt(2)\n");
            flint_abort();
        }

        /* sqrt(2 + sqrt(2)) then squared */
        GR_MUST_SUCCEED(gr_add_ui(e, a, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(e, e, K));
        T = gr_tower_lazy_get_tower(&level, e, K);
        if (level != 2 || gr_tower_degree(T) != 4)
        {
            flint_printf("FAIL: tower of sqrt(2 + sqrt(2))\n");
            gr_tower_print(T);
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_sqr(e, e, K));
        GR_MUST_SUCCEED(gr_sub_ui(e, e, 2, K));
        if (gr_equal(e, a, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(2 + sqrt(2))^2 - 2 != sqrt(2)\n");
            flint_abort();
        }
        gr_tower_lazy_get_tower(&level, e, K);
        if (level != 1)
        {
            flint_printf("FAIL: expected shrink to level 1, got %wd\n", level);
            flint_abort();
        }

        GR_TMP_CLEAR5(a, b, c, d, e, K);
        gr_ctx_clear(K);
    }

    /* Random expressions compared against qqbar arithmetic */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, w;
        gr_ptr a, b, c;
        int flags = n_randint(state, 2) ? GR_TOWER_MERGE_EXPRESS : 0;

        gr_ctx_init_tower_lazy(K, QQ, flags);
        qqbar_init(x); qqbar_init(y); qqbar_init(z); qqbar_init(w);
        GR_TMP_INIT3(a, b, c, K);

        qqbar_randtest(x, state, 1 + n_randint(state, 3), 6);
        qqbar_randtest(y, state, 1 + n_randint(state, 3), 6);

        GR_MUST_SUCCEED(gr_set_other(a, x, QQbar, K));
        GR_MUST_SUCCEED(gr_set_other(b, y, QQbar, K));

        /* c = (a + 1)^2 * b - a / (b^2 + 1) */
        GR_MUST_SUCCEED(gr_add_ui(c, a, 1, K));
        GR_MUST_SUCCEED(gr_sqr(c, c, K));
        GR_MUST_SUCCEED(gr_mul(c, c, b, K));
        {
            gr_ptr t;
            GR_TMP_INIT(t, K);
            GR_MUST_SUCCEED(gr_sqr(t, b, K));
            GR_MUST_SUCCEED(gr_add_ui(t, t, 1, K));
            if (gr_is_zero(t, K) == T_TRUE)   /* y = +/- i */
                GR_MUST_SUCCEED(gr_one(t, K));
            GR_MUST_SUCCEED(gr_div(t, a, t, K));
            GR_MUST_SUCCEED(gr_sub(c, c, t, K));
            GR_TMP_CLEAR(t, K);
        }

        qqbar_add_ui(z, x, 1);
        qqbar_pow_ui(z, z, 2);
        qqbar_mul(z, z, y);
        qqbar_pow_ui(w, y, 2);
        qqbar_add_ui(w, w, 1);
        if (qqbar_is_zero(w))
            qqbar_one(w);
        qqbar_div(w, x, w);
        qqbar_sub(z, z, w);

        GR_MUST_SUCCEED(gr_tower_lazy_get_qqbar(w, c, K));
        if (!qqbar_equal(z, w))
        {
            flint_printf("FAIL: lazy arithmetic vs qqbar\n");
            qqbar_print(z); flint_printf("\n");
            qqbar_print(w); flint_printf("\n");
            flint_abort();
        }

        /* and equality with an independently constructed copy of the same number */
        GR_MUST_SUCCEED(gr_set_other(a, z, QQbar, K));
        if (gr_equal(a, c, K) != T_TRUE)
        {
            flint_printf("FAIL: equality across towers\n");
            gr_println(a, K); gr_println(c, K);
            flint_abort();
        }

        GR_TMP_CLEAR3(a, b, c, K);
        qqbar_clear(x); qqbar_clear(y); qqbar_clear(z); qqbar_clear(w);
        gr_ctx_clear(K);
    }

    /* Transcendental functions: identities which hold in the field of
       rational functions of the generators, numerical separation, and
       identities requiring the relation search. */
    {
        gr_ctx_t K;
        gr_ptr e, p, s, t, u;
        acb_t z, w;

        acb_init(z);
        acb_init(w);
        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_TMP_INIT5(e, p, s, t, u, K);

        GR_MUST_SUCCEED(gr_one(e, K));
        GR_MUST_SUCCEED(gr_exp(e, e, K));           /* e */
        GR_MUST_SUCCEED(gr_pi(p, K));               /* pi */
        GR_MUST_SUCCEED(gr_set_ui(s, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(s, s, K));          /* sqrt(2) */

        /* (pi + e)^2 - pi^2 - 2 pi e == e^2 */
        GR_MUST_SUCCEED(gr_add(t, p, e, K));
        GR_MUST_SUCCEED(gr_sqr(t, t, K));
        GR_MUST_SUCCEED(gr_sqr(u, p, K));
        GR_MUST_SUCCEED(gr_sub(t, t, u, K));
        GR_MUST_SUCCEED(gr_mul(u, p, e, K));
        GR_MUST_SUCCEED(gr_mul_ui(u, u, 2, K));
        GR_MUST_SUCCEED(gr_sub(t, t, u, K));
        GR_MUST_SUCCEED(gr_sqr(u, e, K));
        if (gr_equal(t, u, K) != T_TRUE)
        {
            flint_printf("FAIL: (pi + e)^2\n");
            gr_println(t, K); gr_println(u, K);
            flint_abort();
        }

        /* exp(sqrt(2)) exp(-sqrt(2)) = 1: an identity requiring a relation
           between two generators */
        GR_MUST_SUCCEED(gr_exp(t, s, K));
        GR_MUST_SUCCEED(gr_neg(u, s, K));
        GR_MUST_SUCCEED(gr_exp(u, u, K));
        GR_MUST_SUCCEED(gr_mul(t, t, u, K));
        GR_MUST_SUCCEED(gr_sub_ui(t, t, 1, K));
        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: exp(s2) exp(-s2) - 1\n");
            flint_abort();
        }

        /* the same argument gives the same generator */
        GR_MUST_SUCCEED(gr_exp(t, s, K));
        GR_MUST_SUCCEED(gr_exp(u, s, K));
        if (gr_equal(t, u, K) != T_TRUE)
        {
            flint_printf("FAIL: exp(s2) == exp(s2)\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_one(t, K));
        GR_MUST_SUCCEED(gr_exp(t, t, K));
        if (gr_equal(t, e, K) != T_TRUE)
        {
            flint_printf("FAIL: exp(1) == exp(1)\n");
            flint_abort();
        }

        /* exp(0) = 1, log(1) = 0, log(0) fails */
        GR_MUST_SUCCEED(gr_zero(t, K));
        GR_MUST_SUCCEED(gr_exp(t, t, K));
        if (gr_is_one(t, K) != T_TRUE)
        {
            flint_printf("FAIL: exp(0)\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_one(t, K));
        GR_MUST_SUCCEED(gr_log(t, t, K));
        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: log(1)\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_zero(t, K));
        if (gr_log(t, t, K) != GR_DOMAIN)
        {
            flint_printf("FAIL: log(0)\n");
            flint_abort();
        }

        /* numerical separation: e^2 != exp(3), and division by pi - e */
        GR_MUST_SUCCEED(gr_sqr(t, e, K));
        GR_MUST_SUCCEED(gr_set_ui(u, 3, K));
        GR_MUST_SUCCEED(gr_exp(u, u, K));
        if (gr_equal(t, u, K) != T_FALSE)
        {
            flint_printf("FAIL: e^2 != exp(3)\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_sub(t, p, e, K));
        GR_MUST_SUCCEED(gr_inv(t, t, K));
        GR_MUST_SUCCEED(gr_tower_lazy_get_acb(z, t, 64, K));
        acb_const_pi(w, 64);
        {
            acb_t ee;
            acb_init(ee);
            acb_one(ee); acb_exp(ee, ee, 64);
            acb_sub(w, w, ee, 64);
            acb_inv(w, w, 64);
            acb_clear(ee);
        }
        if (!acb_overlaps(z, w))
        {
            flint_printf("FAIL: 1/(pi - e)\n");
            flint_abort();
        }

        /* log(e + 1) + log(pi) - log(pi (e + 1)) = 0 */
        GR_MUST_SUCCEED(gr_add_ui(t, e, 1, K));
        GR_MUST_SUCCEED(gr_log(t, t, K));
        GR_MUST_SUCCEED(gr_log(u, p, K));
        GR_MUST_SUCCEED(gr_add(t, t, u, K));
        GR_MUST_SUCCEED(gr_add_ui(u, e, 1, K));
        GR_MUST_SUCCEED(gr_mul(u, u, p, K));
        GR_MUST_SUCCEED(gr_log(u, u, K));
        GR_MUST_SUCCEED(gr_sub(t, t, u, K));
        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: log(e+1) + log(pi) - log(pi(e+1))\n");
            flint_abort();
        }

        /* log(e + e^2) - 1 - log(1 + e) = 0 */
        GR_MUST_SUCCEED(gr_sqr(t, e, K));
        GR_MUST_SUCCEED(gr_add(t, t, e, K));
        GR_MUST_SUCCEED(gr_log(t, t, K));
        GR_MUST_SUCCEED(gr_sub_ui(t, t, 1, K));
        GR_MUST_SUCCEED(gr_add_ui(u, e, 1, K));
        GR_MUST_SUCCEED(gr_log(u, u, K));
        GR_MUST_SUCCEED(gr_sub(t, t, u, K));
        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: log(e + e^2) - 1 - log(1 + e)\n");
            flint_abort();
        }

        /* exp(1/2) exp(1/3) = exp(5/6) */
        {
            fmpq_t q;
            fmpq_init(q);
            fmpq_set_si(q, 1, 2);
            GR_MUST_SUCCEED(gr_set_fmpq(t, q, K));
            GR_MUST_SUCCEED(gr_exp(t, t, K));
            fmpq_set_si(q, 1, 3);
            GR_MUST_SUCCEED(gr_set_fmpq(u, q, K));
            GR_MUST_SUCCEED(gr_exp(u, u, K));
            GR_MUST_SUCCEED(gr_mul(t, t, u, K));
            fmpq_set_si(q, 5, 6);
            GR_MUST_SUCCEED(gr_set_fmpq(u, q, K));
            GR_MUST_SUCCEED(gr_exp(u, u, K));
            if (gr_equal(t, u, K) != T_TRUE)
            {
                flint_printf("FAIL: exp(1/2) exp(1/3) = exp(5/6)\n");
                flint_abort();
            }
            fmpq_clear(q);
        }

        /* 3^(log 2 / log 3) = exp((log 2 / log 3) log 3) = 2 */
        GR_MUST_SUCCEED(gr_set_ui(t, 2, K));
        GR_MUST_SUCCEED(gr_log(t, t, K));
        GR_MUST_SUCCEED(gr_set_ui(u, 3, K));
        GR_MUST_SUCCEED(gr_log(u, u, K));
        GR_MUST_SUCCEED(gr_div(t, t, u, K));
        GR_MUST_SUCCEED(gr_mul(t, t, u, K));
        GR_MUST_SUCCEED(gr_exp(t, t, K));
        GR_MUST_SUCCEED(gr_sub_ui(t, t, 2, K));
        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: 3^(log2/log3) = 2\n");
            flint_abort();
        }

        /* exp(pi i) + 1 = 0 and log(-1) = pi i */
        {
            gr_ptr im;
            GR_TMP_INIT(im, K);
            GR_MUST_SUCCEED(gr_set_si(im, -1, K));
            GR_MUST_SUCCEED(gr_sqrt(im, im, K));
            GR_MUST_SUCCEED(gr_mul(t, p, im, K));
            GR_MUST_SUCCEED(gr_exp(u, t, K));
            GR_MUST_SUCCEED(gr_add_ui(u, u, 1, K));
            if (gr_is_zero(u, K) != T_TRUE)
            {
                flint_printf("FAIL: exp(pi i) + 1 = 0\n");
                flint_abort();
            }
            GR_MUST_SUCCEED(gr_set_si(u, -1, K));
            GR_MUST_SUCCEED(gr_log(u, u, K));
            if (gr_equal(u, t, K) != T_TRUE)
            {
                flint_printf("FAIL: log(-1) = pi i\n");
                flint_abort();
            }
            /* and a non-identity of the same shape */
            GR_MUST_SUCCEED(gr_mul_ui(t, t, 2, K));
            if (gr_equal(u, t, K) != T_FALSE)
            {
                flint_printf("FAIL: log(-1) != 2 pi i\n");
                flint_abort();
            }
            GR_TMP_CLEAR(im, K);
        }

        GR_TMP_CLEAR5(e, p, s, t, u, K);
        gr_ctx_clear(K);
        acb_clear(z);
        acb_clear(w);
    }

    /* Polynomial roots requiring new generators: x^2 - e and x^3 - 2 */
    {
        gr_ctx_t K;
        gr_poly_t f;
        gr_vec_t roots;
        fmpz_vec_t mult;
        gr_ptr e, t, u;

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_poly_init(f, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);
        GR_TMP_INIT3(e, t, u, K);

        GR_MUST_SUCCEED(gr_one(e, K));
        GR_MUST_SUCCEED(gr_exp(e, e, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 2, 1, K));
        GR_MUST_SUCCEED(gr_neg(t, e, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, 0, t, K));
        GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));
        if (roots->length != 2)
        {
            flint_printf("FAIL: roots of x^2 - e\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_add(t, gr_vec_entry_ptr(roots, 0, K), gr_vec_entry_ptr(roots, 1, K), K));
        GR_MUST_SUCCEED(gr_sqr(u, gr_vec_entry_ptr(roots, 0, K), K));
        if (gr_is_zero(t, K) != T_TRUE || gr_equal(u, e, K) != T_TRUE)
        {
            flint_printf("FAIL: roots of x^2 - e (values)\n");
            gr_println(gr_vec_entry_ptr(roots, 0, K), K);
            gr_println(gr_vec_entry_ptr(roots, 1, K), K);
            flint_abort();
        }
        /* sqrt(e) computed independently is one of them */
        GR_MUST_SUCCEED(gr_sqrt(t, e, K));
        if (gr_equal(t, gr_vec_entry_ptr(roots, 0, K), K) == T_FALSE && gr_equal(t, gr_vec_entry_ptr(roots, 1, K), K) == T_FALSE)
        {
            flint_printf("FAIL: sqrt(e) is not among the roots of x^2 - e\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_poly_zero(f, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 3, 1, K));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 0, -2, K));
        GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K));
        if (roots->length != 3)
        {
            flint_printf("FAIL: roots of x^3 - 2\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_mul(t, gr_vec_entry_ptr(roots, 0, K), gr_vec_entry_ptr(roots, 1, K), K));
        GR_MUST_SUCCEED(gr_mul(t, t, gr_vec_entry_ptr(roots, 2, K), K));
        GR_MUST_SUCCEED(gr_sub_ui(t, t, 2, K));
        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: product of the roots of x^3 - 2\n");
            gr_println(t, K);
            flint_abort();
        }

        GR_TMP_CLEAR3(e, t, u, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
        gr_poly_clear(f, K);
        gr_ctx_clear(K);
    }

    /* Polynomial roots: random polynomials with coefficients in the lazy
       field, built as products of known linear factors (possibly with
       repetition); the roots must be recovered with multiplicities. */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t K;
        gr_poly_t f, lin;
        gr_vec_t roots;
        fmpz_vec_t mult;
        gr_ptr r, s2, e;
        slong i, n = 1 + n_randint(state, 3), j;
        int status;

        gr_ctx_init_tower_lazy(K, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_poly_init(f, K);
        gr_poly_init(lin, K);
        gr_vec_init(roots, 0, K);
        fmpz_vec_init(mult, 0);
        GR_TMP_INIT3(r, s2, e, K);

        GR_MUST_SUCCEED(gr_set_ui(s2, 2 + n_randint(state, 3), K));
        GR_MUST_SUCCEED(gr_sqrt(s2, s2, K));
        GR_MUST_SUCCEED(gr_set_ui(e, 1 + n_randint(state, 2), K));
        GR_MUST_SUCCEED(gr_exp(e, e, K));

        GR_MUST_SUCCEED(gr_poly_one(f, K));
        for (i = 0; i < n; i++)
        {
            /* root: a small rational, a multiple of sqrt(2) or of e */
            switch (n_randint(state, 3))
            {
                case 0: GR_MUST_SUCCEED(gr_set_si(r, (slong) n_randint(state, 5) - 2, K)); break;
                case 1: GR_MUST_SUCCEED(gr_mul_si(r, s2, (slong) n_randint(state, 3) - 1, K)); break;
                default: GR_MUST_SUCCEED(gr_mul_si(r, e, (slong) n_randint(state, 3) - 1, K)); break;
            }
            GR_MUST_SUCCEED(gr_poly_zero(lin, K));
            GR_MUST_SUCCEED(gr_poly_set_coeff_si(lin, 1, 1, K));
            GR_MUST_SUCCEED(gr_neg(r, r, K));
            GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(lin, 0, r, K));
            GR_MUST_SUCCEED(gr_poly_mul(f, f, lin, K));
            if (n_randint(state, 3) == 0)
                GR_MUST_SUCCEED(gr_poly_mul(f, f, lin, K));   /* repeated root */
        }
        /* a non-monic leading coefficient */
        GR_MUST_SUCCEED(gr_poly_mul_scalar(f, f, s2, K));

        status = gr_poly_roots(roots, mult, f, 0, K);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: roots status %d\n", status);
            gr_poly_print(f, K); flint_printf("\n");
            flint_abort();
        }

        /* the roots are roots, and their multiplicities add up to the degree */
        {
            slong total = 0;
            for (i = 0; i < roots->length; i++)
            {
                GR_MUST_SUCCEED(gr_poly_evaluate(r, f, gr_vec_entry_ptr(roots, i, K), K));
                if (gr_is_zero(r, K) != T_TRUE)
                {
                    flint_printf("FAIL: root %wd is not a root\n", i);
                    gr_poly_print(f, K); flint_printf("\n");
                    gr_println(gr_vec_entry_ptr(roots, i, K), K);
                    flint_abort();
                }
                total += fmpz_get_si(mult->entries + i);
                for (j = 0; j < i; j++)
                {
                    if (gr_equal(gr_vec_entry_ptr(roots, i, K), gr_vec_entry_ptr(roots, j, K), K) != T_FALSE)
                    {
                        flint_printf("FAIL: roots %wd and %wd not distinct\n", i, j);
                        flint_abort();
                    }
                }
            }
            if (total != f->length - 1)
            {
                flint_printf("FAIL: multiplicities %wd != degree %wd\n", total, f->length - 1);
                gr_poly_print(f, K); flint_printf("\n");
                flint_abort();
            }
        }

        GR_TMP_CLEAR3(r, s2, e, K);
        gr_vec_clear(roots, K);
        fmpz_vec_clear(mult);
        gr_poly_clear(f, K);
        gr_poly_clear(lin, K);
        gr_ctx_clear(K);
    }

    /* Principal roots next to the branch cut: x = (e + i)^n with
       e = sqrt(2) - 1.41421356237309504880168872 (about 4e-27) has an
       enclosure straddling the negative real axis at low precision; the
       n-th root must still be e + i. */
    {
        gr_ptr e, w, x, r;
        ulong n;
        fmpq_t c;

        gr_ctx_init_tower_lazy(K, QQ, 0);
        GR_TMP_INIT4(e, w, x, r, K);
        fmpq_init(c);
        fmpz_set_str(fmpq_numref(c), "141421356237309504880168872", 10);
        fmpz_set_str(fmpq_denref(c), "100000000000000000000000000", 10);
        GR_MUST_SUCCEED(gr_set_ui(e, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(e, e, K));
        GR_MUST_SUCCEED(gr_sub_fmpq(e, e, c, K));
        GR_MUST_SUCCEED(gr_i(w, K));
        GR_MUST_SUCCEED(gr_add(w, w, e, K));

        for (n = 2; n <= 5; n++)
        {
            int status;
            gr_ptr wn, zeta;
            qqbar_t q;

            /* w_n = (e + i) zeta_{4n}^(2-n) has argument pi/n - tiny, so
               that w_n^n straddles the cut and w_n is its principal root */
            GR_TMP_INIT2(wn, zeta, K);
            qqbar_init(q);
            qqbar_root_of_unity(q, (slong) 2 - (slong) n, 4 * n);
            GR_MUST_SUCCEED(gr_set_other(zeta, q, QQbar, K));
            GR_MUST_SUCCEED(gr_mul(wn, w, zeta, K));
            qqbar_clear(q);

            GR_MUST_SUCCEED(gr_pow_ui(x, wn, n, K));
            status = (n == 2) ? gr_sqrt(r, x, K) : gr_tower_lazy_root_ui(r, x, n, K);
            if (status != GR_SUCCESS && status != GR_UNABLE)
            {
                flint_printf("FAIL: root %wu status %d\n", n, status);
                flint_abort();
            }
            if (status == GR_SUCCESS && gr_equal(r, wn, K) != T_TRUE)
            {
                flint_printf("FAIL: non-principal root %wu\n", n);
                gr_println(r, K);
                flint_abort();
            }
            GR_TMP_CLEAR2(wn, zeta, K);
        }

        fmpq_clear(c);
        GR_TMP_CLEAR4(e, w, x, r, K);
        gr_ctx_clear(K);
    }

    /* Transfer between two lazy contexts and string round trips must be
       exact even for roots that are closer than their printed approximations
       would suggest: x^20 - 2 (101 x - 1)^2 has two real roots about
       1.3e-22 apart near 1/101. */
    {
        gr_ctx_t K1, K2;
        gr_poly_t f, lin;
        gr_vec_t roots;
        fmpz_vec_t mult;
        gr_ptr y1, y2, z;
        slong i, j, nreal;
        gr_ptr xs[2];
        char * str;

        gr_ctx_init_tower_lazy(K1, QQ, 0);
        gr_ctx_init_tower_lazy(K2, QQ, 0);

        gr_poly_init(f, K1);
        gr_poly_init(lin, K1);
        gr_vec_init(roots, 0, K1);
        fmpz_vec_init(mult, 0);
        GR_TMP_INIT3(y1, y2, z, K2);

        /* f = x^20 - 2 (101 x - 1)^2 */
        GR_MUST_SUCCEED(gr_poly_zero(lin, K1));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(lin, 1, 101, K1));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(lin, 0, -1, K1));
        GR_MUST_SUCCEED(gr_poly_mul(f, lin, lin, K1));
        GR_MUST_SUCCEED(gr_poly_neg(f, f, K1));
        GR_MUST_SUCCEED(gr_poly_add(f, f, f, K1));
        GR_MUST_SUCCEED(gr_poly_set_coeff_si(f, 20, 1, K1));
        GR_MUST_SUCCEED(gr_poly_roots(roots, mult, f, 0, K1));

        /* pick the two real roots closest to 1/101 */
        nreal = 0;
        for (i = 0; i < roots->length && nreal < 2; i++)
        {
            gr_ptr r = gr_vec_entry_ptr(roots, i, K1);
            gr_ptr t;
            GR_TMP_INIT(t, K1);
            GR_MUST_SUCCEED(gr_mul_si(t, r, 101, K1));
            GR_MUST_SUCCEED(gr_sub_ui(t, t, 1, K1));
            GR_MUST_SUCCEED(gr_abs(t, t, K1));
            GR_MUST_SUCCEED(gr_mul_2exp_si(t, t, 60, K1));
            /* real root with |101 r - 1| < 2^-60 */
            if (gr_tower_lazy_is_real(r, K1) == T_TRUE)
            {
                gr_ptr one;
                GR_TMP_INIT(one, K1);
                GR_MUST_SUCCEED(gr_one(one, K1));
                if (gr_lt(t, one, K1) == T_TRUE)
                    xs[nreal++] = r;
                GR_TMP_CLEAR(one, K1);
            }
            GR_TMP_CLEAR(t, K1);
        }

        if (nreal != 2)
        {
            flint_printf("FAIL: expected two close roots near 1/101, found %wd\n", nreal);
            gr_poly_print(f, K1); flint_printf("\n");
            flint_abort();
        }

        /* transfer both roots into K2: they must stay distinct and map back */
        GR_MUST_SUCCEED(gr_set_other(y1, xs[0], K1, K2));
        GR_MUST_SUCCEED(gr_set_other(y2, xs[1], K1, K2));

        if (gr_equal(y1, y2, K2) != T_FALSE)
        {
            flint_printf("FAIL: close roots became equal after transfer\n");
            gr_println(y1, K2); gr_println(y2, K2);
            flint_abort();
        }

        for (j = 0; j < 2; j++)
        {
            gr_ptr back;
            GR_TMP_INIT(back, K1);
            GR_MUST_SUCCEED(gr_set_other(back, j == 0 ? y1 : y2, K2, K1));
            if (gr_equal(back, xs[j], K1) != T_TRUE)
            {
                flint_printf("FAIL: transfer K1 -> K2 -> K1 is not the identity (root %wd)\n", j);
                gr_println(xs[j], K1); gr_println(back, K1);
                flint_abort();
            }
            GR_TMP_CLEAR(back, K1);

            /* string round trip within K1 */
            GR_MUST_SUCCEED(gr_get_str(&str, xs[j], K1));
            GR_MUST_SUCCEED(gr_set_str(z, str, K2));
            /* compare via the transferred copies */
            if (gr_equal(z, j == 0 ? y1 : y2, K2) != T_TRUE)
            {
                flint_printf("FAIL: string round trip of root %wd: %s\n", j, str);
                gr_println(z, K2);
                flint_abort();
            }
            flint_free(str);
        }

        GR_TMP_CLEAR3(y1, y2, z, K2);
        gr_vec_clear(roots, K1);
        fmpz_vec_clear(mult);
        gr_poly_clear(f, K1);
        gr_poly_clear(lin, K1);
        gr_ctx_clear(K1);
        gr_ctx_clear(K2);
    }

    gr_ctx_clear(QQ);
    gr_ctx_clear(QQbar);

    TEST_FUNCTION_END(state);
}
