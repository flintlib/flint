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
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "gr_poly.h"
#include "acb.h"
#include "qqbar.h"
#include "gr_vec.h"
#include "ulong_extras.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

/* Q(sqrt(p_1), ..., sqrt(p_n)) for the first n primes */
static void
_multiquadratic(gr_tower_t T, gr_ctx_t QQ, slong n)
{
    ulong p = 2;
    slong k;

    gr_tower_init(T, QQ);

    for (k = 0; k < n; k++)
    {
        gr_ctx_struct * top = gr_tower_field(T);
        gr_ptr x;
        GR_TMP_INIT(x, top);
        GR_MUST_SUCCEED(gr_set_ui(x, p, top));
        GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
        GR_TMP_CLEAR(x, top);
        p = n_nextprime(p, 1);
    }
}

TEST_FUNCTION_START(gr_tower_algebraic, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    /* Multiquadratic fields: ring axioms, evaluation, and the identity
       sqrt(2)*sqrt(3)*...*sqrt(p_n) = sqrt(2*3*...*p_n) via qqbar */
    for (iter = 1; iter <= 4; iter++)
    {
        gr_tower_t T;
        gr_ctx_t K;
        gr_vec_t gens;
        gr_ptr prod;
        qqbar_t q1, q2;
        acb_t z1, z2;
        slong k;
        ulong N = 1, p = 2;

        _multiquadratic(T, QQ, iter);

        if (gr_tower_degree(T) != (WORD(1) << iter))
        {
            flint_printf("FAIL: degree\n");
            flint_abort();
        }

        gr_ctx_init_tower_field(K, T);
        gr_test_ring(K, 5 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);

        gr_vec_init(gens, 0, K);
        GR_MUST_SUCCEED(gr_gens(gens, K));
        GR_TMP_INIT(prod, K);
        GR_MUST_SUCCEED(gr_one(prod, K));
        for (k = 0; k < iter; k++)
        {
            GR_MUST_SUCCEED(gr_mul(prod, prod, gr_vec_entry_ptr(gens, k, K), K));
            N *= p;
            p = n_nextprime(p, 1);
        }

        qqbar_init(q1);
        qqbar_init(q2);
        acb_init(z1);
        acb_init(z2);

        GR_MUST_SUCCEED(gr_tower_get_qqbar(q1, prod, T));
        qqbar_set_ui(q2, N);
        qqbar_sqrt(q2, q2);

        if (!qqbar_equal(q1, q2))
        {
            flint_printf("FAIL: product of square roots\n");
            gr_tower_print(T);
            qqbar_print(q1); flint_printf("\n");
            qqbar_print(q2); flint_printf("\n");
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_tower_get_acb(z1, prod, 200, T));
        qqbar_get_acb(z2, q2, 200);
        if (!acb_overlaps(z1, z2) || acb_rel_accuracy_bits(z1) < 190)
        {
            flint_printf("FAIL: enclosure\n");
            acb_printn(z1, 50, 0); flint_printf("\n");
            acb_printn(z2, 50, 0); flint_printf("\n");
            flint_abort();
        }

        qqbar_clear(q1);
        qqbar_clear(q2);
        acb_clear(z1);
        acb_clear(z2);
        GR_TMP_CLEAR(prod, K);
        gr_vec_clear(gens, K);
        gr_ctx_clear(K);
        gr_tower_clear(T);
    }

    /* Dynamic evaluation: adjoin sqrt(2), then the qqbar sqrt(8), whose
       minimal polynomial x^2 - 8 factors over Q(sqrt(2)). The zero test
       for sqrt(8) - 2 sqrt(2) must discover this. */
    {
        gr_tower_t T;
        gr_ctx_t K;
        gr_vec_t gens;
        gr_ptr t;
        qqbar_t q;

        gr_tower_init(T, QQ);
        {
            gr_ptr x;
            GR_TMP_INIT(x, QQ);
            GR_MUST_SUCCEED(gr_set_ui(x, 2, QQ));
            GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, "s2"));
            GR_TMP_CLEAR(x, QQ);
        }

        qqbar_init(q);
        qqbar_set_ui(q, 8);
        qqbar_sqrt(q, q);
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, q, "s8"));

        if (gr_tower_degree(T) != 4)
        {
            flint_printf("FAIL: expected degree 4 before refinement\n");
            flint_abort();
        }

        gr_ctx_init_tower_field(K, T);
        gr_vec_init(gens, 0, K);
        GR_MUST_SUCCEED(gr_gens(gens, K));
        GR_TMP_INIT(t, K);

        /* t = s8 - 2 s2 */
        GR_MUST_SUCCEED(gr_mul_si(t, gr_vec_entry_ptr(gens, 0, K), 2, K));
        GR_MUST_SUCCEED(gr_sub(t, gr_vec_entry_ptr(gens, 1, K), t, K));

        if (gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: sqrt(8) - 2 sqrt(2) != 0\n");
            gr_tower_print(T);
            gr_println(t, K);
            flint_abort();
        }

        if (gr_tower_degree(T) != 2 || gr_tower_step_degree(T, 2) != 1)
        {
            flint_printf("FAIL: expected refinement to degree 2\n");
            gr_tower_print(T);
            flint_abort();
        }

        /* the stale generator vector should still be usable */
        GR_MUST_SUCCEED(gr_add(t, gr_vec_entry_ptr(gens, 0, K), gr_vec_entry_ptr(gens, 1, K), K));
        GR_MUST_SUCCEED(gr_sqr(t, t, K));
        if (gr_sub_ui(t, t, 18, K) != GR_SUCCESS || gr_is_zero(t, K) != T_TRUE)
        {
            flint_printf("FAIL: (s2 + s8)^2 != 18\n");
            gr_println(t, K);
            flint_abort();
        }

        gr_test_ring(K, 5 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);

        qqbar_clear(q);
        GR_TMP_CLEAR(t, K);
        gr_vec_clear(gens, K);
        gr_ctx_clear(K);
        gr_tower_clear(T);
    }

    /* Copying a tower, and charpoly vs minpoly */
    {
        gr_tower_t T, T2;
        gr_ctx_t K, K2;
        gr_vec_t gens;
        gr_ptr e, f;
        gr_poly_t chi, mp;
        qqbar_t q;

        _multiquadratic(T, QQ, 3);
        gr_tower_init(T2, QQ);
        gr_tower_set(T2, T);

        if (gr_tower_degree(T2) != 8 || gr_tower_length(T2) != 3)
        {
            flint_printf("FAIL: tower copy\n");
            flint_abort();
        }

        gr_ctx_init_tower_field(K, T);
        gr_ctx_init_tower_field(K2, T2);
        gr_vec_init(gens, 0, K);
        GR_MUST_SUCCEED(gr_gens(gens, K));
        GR_TMP_INIT(e, K);
        GR_TMP_INIT(f, K2);

        /* e = s2 + s3 + s5; the same polynomial data is a valid element of the copy */
        GR_MUST_SUCCEED(gr_add(e, gr_vec_entry_ptr(gens, 0, K), gr_vec_entry_ptr(gens, 1, K), K));
        GR_MUST_SUCCEED(gr_add(e, e, gr_vec_entry_ptr(gens, 2, K), K));
        GR_MUST_SUCCEED(gr_set(f, e, K2));
        GR_MUST_SUCCEED(gr_sqr(f, f, K2));

        gr_poly_init(chi, QQ);
        gr_poly_init(mp, QQ);
        GR_MUST_SUCCEED(gr_tower_charpoly(chi, e, T));
        GR_MUST_SUCCEED(_gr_tower_annihilating_poly(mp, e, T));

        /* s2 + s3 + s5 is a primitive element of degree 8 */
        if (gr_poly_length(chi, QQ) != 9 || gr_poly_equal(chi, mp, QQ) != T_TRUE)
        {
            flint_printf("FAIL: charpoly / minpoly\n");
            GR_MUST_SUCCEED(gr_poly_print(chi, QQ)); flint_printf("\n");
            GR_MUST_SUCCEED(gr_poly_print(mp, QQ)); flint_printf("\n");
            flint_abort();
        }

        qqbar_init(q);
        /* (s2+s3+s5)^2 = 10 + 2 s6 + 2 s10 + 2 s15 has degree 4 */
        GR_MUST_SUCCEED(gr_tower_get_qqbar(q, f, T2));
        if (qqbar_degree(q) != 4)
        {
            flint_printf("FAIL: degree of (s2+s3+s5)^2\n");
            flint_abort();
        }
        qqbar_clear(q);

        gr_poly_clear(chi, QQ);
        gr_poly_clear(mp, QQ);
        GR_TMP_CLEAR(e, K);
        GR_TMP_CLEAR(f, K2);
        gr_vec_clear(gens, K);
        gr_ctx_clear(K);
        gr_ctx_clear(K2);
        gr_tower_clear(T);
        gr_tower_clear(T2);
    }

    /* Random qqbar round trips: adjoin a random algebraic number, and
       check that get_qqbar of random polynomial expressions agrees with
       qqbar arithmetic. */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        gr_ctx_t K;
        qqbar_t x, y, z, w;
        gr_vec_t gens;
        gr_ptr e, f, g;
        slong deg;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(w);

        deg = 1 + n_randint(state, 4);
        qqbar_randtest(x, state, deg, 10);
        qqbar_randtest(y, state, 2, 10);

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, x, "x"));
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, y, "y"));

        gr_ctx_init_tower_field(K, T);
        gr_vec_init(gens, 0, K);
        GR_MUST_SUCCEED(gr_gens(gens, K));
        GR_TMP_INIT3(e, f, g, K);

        /* e = x*y + 3x - 2, compare with qqbar */
        GR_MUST_SUCCEED(gr_mul(e, gr_vec_entry_ptr(gens, 0, K), gr_vec_entry_ptr(gens, 1, K), K));
        GR_MUST_SUCCEED(gr_mul_si(f, gr_vec_entry_ptr(gens, 0, K), 3, K));
        GR_MUST_SUCCEED(gr_add(e, e, f, K));
        GR_MUST_SUCCEED(gr_sub_ui(e, e, 2, K));

        qqbar_mul(z, x, y);
        qqbar_mul_si(w, x, 3);
        qqbar_add(z, z, w);
        qqbar_sub_ui(z, z, 2);

        GR_MUST_SUCCEED(gr_tower_get_qqbar(w, e, T));

        if (!qqbar_equal(z, w))
        {
            flint_printf("FAIL: qqbar round trip\n");
            gr_tower_print(T);
            qqbar_print(z); flint_printf("\n");
            qqbar_print(w); flint_printf("\n");
            flint_abort();
        }

        /* Division: e / (e^2 + 1) * (e^2 + 1) == e */
        GR_MUST_SUCCEED(gr_sqr(f, e, K));
        GR_MUST_SUCCEED(gr_add_ui(f, f, 1, K));
        GR_MUST_SUCCEED(gr_div(g, e, f, K));
        GR_MUST_SUCCEED(gr_mul(g, g, f, K));
        if (gr_equal(g, e, K) != T_TRUE)
        {
            flint_printf("FAIL: division\n");
            gr_tower_print(T);
            flint_abort();
        }

        /* Zero test of a nonzero element should be decided */
        if (gr_is_zero(f, K) != T_FALSE)
        {
            flint_printf("FAIL: e^2 + 1 == 0 ?\n");
            gr_tower_print(T);
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(w);
        GR_TMP_CLEAR3(e, f, g, K);
        gr_vec_clear(gens, K);
        gr_ctx_clear(K);
        gr_tower_clear(T);
    }

    /* Redundant generators: adjoin x, then y = p(x) as a qqbar. The
       minimal polynomial of y over Q(x) has a linear factor, which must be
       found when testing y - p(x) == 0; the tower then has degree deg(x). */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        gr_ctx_t K;
        qqbar_t x, y;
        fmpq_poly_t pol;
        gr_vec_t gens;
        gr_ptr t, u;
        slong deg;

        qqbar_init(x);
        qqbar_init(y);
        fmpq_poly_init(pol);

        deg = 2 + n_randint(state, 4);
        do {
            qqbar_randtest(x, state, deg, 8);
        } while (qqbar_degree(x) < 2);

        fmpq_poly_randtest(pol, state, 1 + n_randint(state, deg), 6);
        qqbar_evaluate_fmpq_poly(y, pol, x);

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, x, "x"));
        GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, y, "y"));

        gr_ctx_init_tower_field(K, T);
        gr_vec_init(gens, 0, K);
        GR_MUST_SUCCEED(gr_gens(gens, K));
        GR_TMP_INIT2(t, u, K);

        /* t = p(x) */
        {
            gr_poly_t pg;
            gr_poly_init(pg, K);
            GR_MUST_SUCCEED(gr_poly_set_fmpq_poly(pg, pol, K));
            GR_MUST_SUCCEED(gr_poly_evaluate(t, pg, gr_vec_entry_ptr(gens, 0, K), K));
            gr_poly_clear(pg, K);
        }

        if (gr_equal(t, gr_vec_entry_ptr(gens, 1, K), K) != T_TRUE)
        {
            flint_printf("FAIL: y != p(x)\n");
            gr_tower_print(T);
            flint_abort();
        }

        if (gr_tower_degree(T) != qqbar_degree(x))
        {
            flint_printf("FAIL: expected degree %wd after refinement, got %wd\n", qqbar_degree(x), gr_tower_degree(T));
            gr_tower_print(T);
            flint_abort();
        }

        /* the ring must still behave after refinement */
        gr_test_ring(K, 3 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);

        /* and expressions in the (stale) generators still evaluate correctly */
        GR_MUST_SUCCEED(gr_mul(u, gr_vec_entry_ptr(gens, 1, K), gr_vec_entry_ptr(gens, 0, K), K));
        GR_MUST_SUCCEED(gr_mul(t, t, gr_vec_entry_ptr(gens, 0, K), K));
        if (gr_equal(t, u, K) != T_TRUE)
        {
            flint_printf("FAIL: x*y != x*p(x)\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        fmpq_poly_clear(pol);
        GR_TMP_CLEAR2(t, u, K);
        gr_vec_clear(gens, K);
        gr_ctx_clear(K);
        gr_tower_clear(T);
    }

    /* Random towers of random algebraic numbers under the generic ring tests,
       which exercise inversion (and hence refinement) on random elements. */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        gr_ctx_t K;
        qqbar_t x;
        slong k, n = 1 + n_randint(state, 2);

        qqbar_init(x);
        gr_tower_init(T, QQ);

        for (k = 0; k < n; k++)
        {
            if (n_randint(state, 2))
            {
                qqbar_randtest(x, state, 1 + n_randint(state, 3), 6);
                GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, x, NULL));
            }
            else
            {
                gr_ctx_struct * top = gr_tower_field(T);
                gr_ptr e;
                GR_TMP_INIT(e, top);
                GR_MUST_SUCCEED(gr_randtest(e, state, top));
                if (gr_tower_adjoin_root_ui(T, e, 2 + n_randint(state, 2), NULL) != GR_SUCCESS)
                {
                    /* e may be zero or the enclosure may fail to certify; skip */
                }
                GR_TMP_CLEAR(e, top);
            }
        }

        gr_ctx_init_tower_field(K, T);
        gr_test_ring(K, 3 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);
        gr_ctx_clear(K);

        qqbar_clear(x);
        gr_tower_clear(T);
    }

    /* Square roots in multiquadratic towers, without lattice reduction */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        gr_ctx_struct * top;
        gr_ptr x, y, r;
        acb_t zx, zr;
        int status;

        _multiquadratic(T, QQ, 1 + n_randint(state, 4));
        top = gr_tower_field(T);
        GR_TMP_INIT3(x, y, r, top);
        acb_init(zx);
        acb_init(zr);

        GR_MUST_SUCCEED(gr_randtest(y, state, top));
        GR_MUST_SUCCEED(gr_sqr(x, y, top));

        /* sqrt(y^2) is the principal root: +-y, with the right sign */
        status = gr_tower_sqrt(r, x, T);
        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: sqrt of a square (status %d)\n", status);
            gr_tower_print(T);
            gr_println(y, top);
            flint_abort();
        }
        if (gr_tower_equal(r, y, T) != T_TRUE)
        {
            GR_MUST_SUCCEED(gr_neg(r, r, top));
            if (gr_tower_equal(r, y, T) != T_TRUE)
            {
                flint_printf("FAIL: sqrt(y^2) != +-y\n");
                flint_abort();
            }
            GR_MUST_SUCCEED(gr_neg(r, r, top));
        }
        GR_MUST_SUCCEED(gr_tower_get_acb(zx, x, 64, T));
        GR_MUST_SUCCEED(gr_tower_get_acb(zr, r, 64, T));
        acb_sqrt(zx, zx, 64);
        if (!acb_overlaps(zx, zr))
        {
            flint_printf("FAIL: sqrt sign\n");
            flint_abort();
        }

        /* a random element: if a root is found, it is a root */
        GR_MUST_SUCCEED(gr_add_ui(x, x, 1 + n_randint(state, 3), top));
        status = gr_tower_sqrt(r, x, T);
        if (status == GR_SUCCESS)
        {
            GR_MUST_SUCCEED(gr_sqr(r, r, top));
            if (gr_tower_equal(r, x, T) != T_TRUE)
            {
                flint_printf("FAIL: sqrt(x)^2 != x\n");
                flint_abort();
            }
        }
        else if (status != GR_DOMAIN)
        {
            flint_printf("FAIL: sqrt status %d\n", status);
            flint_abort();
        }

        GR_TMP_CLEAR3(x, y, r, top);
        acb_clear(zx);
        acb_clear(zr);
        gr_tower_clear(T);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
