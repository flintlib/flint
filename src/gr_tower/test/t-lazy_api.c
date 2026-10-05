/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Coverage of the lazy field methods which the other tests do not
   exercise: inverse trigonometric and hyperbolic functions at complex
   points against acb, arg and sgn, comparisons and realness against
   balls, conversions to doubles, beta and Hurwitz zeta values, and
   conversions between contexts (polynomial factorization is covered by
   t-factor.c). */

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

#define CHECK(cond, msg) do { if (!(cond)) { flint_printf("FAIL: %s\n", msg); flint_abort(); } } while (0)

#define PREC 128

/* a random algebraic number with a few generators: q + sqrt(a) + c i */
static void
_rand_alg(gr_ptr x, acb_t z, flint_rand_t state, gr_ctx_t K)
{
    fmpq_t q;
    gr_ptr t;
    slong a = 2 + n_randint(state, 7), c = (slong) n_randint(state, 5) - 2;
    acb_t w;

    fmpq_init(q);
    acb_init(w);
    GR_TMP_INIT(t, K);

    fmpq_set_si(q, (slong) n_randint(state, 21) - 10, 1 + n_randint(state, 4));
    GR_MUST_SUCCEED(gr_set_fmpq(x, q, K));
    acb_set_fmpq(z, q, PREC);

    if (n_randint(state, 2))
    {
        GR_MUST_SUCCEED(gr_set_si(t, a, K));
        GR_MUST_SUCCEED(gr_sqrt(t, t, K));
        GR_MUST_SUCCEED(gr_add(x, x, t, K));
        acb_set_si(w, a);
        acb_sqrt(w, w, PREC);
        acb_add(z, z, w, PREC);
    }

    if (c != 0 && n_randint(state, 2))
    {
        GR_MUST_SUCCEED(gr_i(t, K));
        GR_MUST_SUCCEED(gr_mul_si(t, t, c, K));
        GR_MUST_SUCCEED(gr_add(x, x, t, K));
        acb_onei(w);
        acb_mul_si(w, w, c, PREC);
        acb_add(z, z, w, PREC);
    }

    fmpq_clear(q);
    acb_clear(w);
    GR_TMP_CLEAR(t, K);
}

static int
_overlaps(gr_srcptr x, const acb_t v, gr_ctx_t K)
{
    acb_t w;
    int ok;
    acb_init(w);
    GR_MUST_SUCCEED(gr_tower_lazy_get_acb(w, x, PREC, K));
    ok = acb_overlaps(w, v);
    acb_clear(w);
    return ok;
}

TEST_FUNCTION_START(gr_tower_lazy_api, state)
{
    gr_ctx_t QQ, K, R;
    slong iter;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_tower_lazy(K, QQ, 0);
    gr_ctx_init_tower_lazy_view(R, K, GR_TOWER_LAZY_REAL);

    /* inverse functions at complex algebraic points (including points
       on the branch cuts), arg and sgn, against acb */
    {
    slong successes = 0, total = 60 * flint_test_multiplier();
    for (iter = 0; iter < total; iter++)
    {
        gr_ptr x, y;
        acb_t z, v;
        int which = n_randint(state, 8), status;

        GR_TMP_INIT2(x, y, K);
        acb_init(z);
        acb_init(v);
        _rand_alg(x, z, state, K);

        switch (which)
        {
            case 0: status = gr_asin(y, x, K); acb_asin(v, z, PREC); break;
            case 1: status = gr_acos(y, x, K); acb_acos(v, z, PREC); break;
            case 2: status = gr_atan(y, x, K); acb_atan(v, z, PREC); break;
            case 3: status = gr_asinh(y, x, K); acb_asinh(v, z, PREC); break;
            case 4: status = gr_acosh(y, x, K); acb_acosh(v, z, PREC); break;
            case 5: status = gr_atanh(y, x, K); acb_atanh(v, z, PREC); break;
            case 6: status = gr_arg(y, x, K); acb_arg(acb_realref(v), z, PREC); arb_zero(acb_imagref(v)); break;
            default: status = gr_sgn(y, x, K); acb_sgn(v, z, PREC); break;
        }

        if (status == GR_SUCCESS)
        {
            successes++;
            if (!_overlaps(y, v, K))
            {
                flint_printf("FAIL: function %d\n", which);
                flint_printf("x = "); gr_println(x, K);
                flint_printf("y = "); gr_println(y, K);
                acb_printn(v, 30, 0); flint_printf("\n");
                flint_abort();
            }
        }
        else if (status == GR_DOMAIN)
        {
            CHECK(!acb_is_finite(v), "domain error at a finite value");
        }

        GR_TMP_CLEAR2(x, y, K);
        acb_clear(z);
        acb_clear(v);
    }
    CHECK(total < 8 || successes >= total / 2, "most inverse function values are computed");
    }

    /* comparisons, signs and realness against balls */
    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        gr_ptr x, y;
        acb_t z, w;
        int c, c2, status;
        truth_t real;

        GR_TMP_INIT2(x, y, K);
        acb_init(z);
        acb_init(w);
        _rand_alg(x, z, state, K);
        _rand_alg(y, w, state, K);

        real = gr_tower_lazy_is_real(x, K);
        CHECK(real != T_UNKNOWN, "realness decided");
        CHECK((real == T_TRUE) == arb_is_zero(acb_imagref(z)), "realness");

        status = gr_cmpabs(&c, x, y, K);
        if (status == GR_SUCCESS)
        {
            arb_t ax, ay;
            arb_init(ax);
            arb_init(ay);
            acb_abs(ax, z, PREC);
            acb_abs(ay, w, PREC);
            if (c != 0)
                CHECK((c < 0) == (arb_lt(ax, ay) != 0) || arb_overlaps(ax, ay), "cmpabs");
            else
                CHECK(arb_overlaps(ax, ay), "cmpabs equal");
            arb_clear(ax);
            arb_clear(ay);
        }

        if (real == T_TRUE && gr_tower_lazy_is_real(y, K) == T_TRUE)
        {
            status = gr_cmp(&c, x, y, K);
            CHECK(status == GR_SUCCESS, "cmp of real elements");
            status = gr_tower_lazy_cmp(&c2, x, y, K);
            CHECK(status == GR_SUCCESS && c2 == c, "lazy_cmp");
            if (c != 0)
                CHECK((c < 0) == (arb_lt(acb_realref(z), acb_realref(w)) != 0) || arb_overlaps(acb_realref(z), acb_realref(w)), "cmp");
            else
                CHECK(arb_overlaps(acb_realref(z), acb_realref(w)), "cmp equal");

            /* the sign in the real view */
            {
                gr_ptr s;
                GR_TMP_INIT(s, R);
                GR_MUST_SUCCEED(gr_set_other(s, x, K, R));
                GR_MUST_SUCCEED(gr_tower_lazy_sgn(s, s, R));
                CHECK(gr_is_zero(s, R) == (arb_is_zero(acb_realref(z)) ? T_TRUE : T_FALSE), "sgn zero");
                if (arb_is_positive(acb_realref(z)))
                    CHECK(gr_is_one(s, R) == T_TRUE, "sgn positive");
                if (arb_is_negative(acb_realref(z)))
                    CHECK(gr_is_neg_one(s, R) == T_TRUE, "sgn negative");
                GR_TMP_CLEAR(s, R);
            }
        }
        else
        {
            CHECK(gr_set_other(y, x, K, R) == GR_DOMAIN || real == T_TRUE, "nonreal element in a real view");
        }

        /* conversion to a double */
        {
            double d;
            status = gr_tower_lazy_get_d(&d, x, K);
            if (real == T_TRUE)
            {
                arb_t t;
                arb_init(t);
                CHECK(status == GR_SUCCESS, "get_d of a real element");
                arb_set_d(t, d);
                arb_add_error_2exp_si(t, -40);
                CHECK(arb_overlaps(t, acb_realref(z)), "get_d value");
                arb_clear(t);
            }
            else
                CHECK(status == GR_DOMAIN, "get_d of a nonreal element");
        }

        GR_TMP_CLEAR2(x, y, K);
        acb_clear(z);
        acb_clear(w);
    }

    /* beta and Hurwitz zeta values: B(1/2, 1/2) = pi, B(2, 3) = 1/12,
       zeta(2, 1) = pi^2/6, zeta(2, 1/2) = pi^2/2, zeta(s, a) numerically */
    {
        gr_ptr x, y, z, pi;
        fmpq_t q;
        acb_t s, a, v;
        int status;

        GR_TMP_INIT4(x, y, z, pi, K);
        fmpq_init(q);
        acb_init(s);
        acb_init(a);
        acb_init(v);
        GR_MUST_SUCCEED(gr_pi(pi, K));

        fmpq_set_si(q, 1, 2);
        GR_MUST_SUCCEED(gr_set_fmpq(x, q, K));
        status = gr_tower_lazy_beta(z, x, x, K);
        CHECK(status == GR_SUCCESS && gr_equal(z, pi, K) == T_TRUE, "B(1/2, 1/2) = pi");

        GR_MUST_SUCCEED(gr_set_ui(x, 2, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 3, K));
        status = gr_tower_lazy_beta(z, x, y, K);
        fmpq_set_si(q, 1, 12);
        GR_MUST_SUCCEED(gr_set_fmpq(y, q, K));
        CHECK(status == GR_SUCCESS && gr_equal(z, y, K) == T_TRUE, "B(2, 3) = 1/12");

        GR_MUST_SUCCEED(gr_set_ui(x, 2, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 1, K));
        status = gr_tower_lazy_hurwitz_zeta(z, x, y, K);
        GR_MUST_SUCCEED(gr_mul(y, pi, pi, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 6, K));
        CHECK(status == GR_SUCCESS && gr_equal(z, y, K) == T_TRUE, "zeta(2, 1) = pi^2/6");

        fmpq_set_si(q, 1, 2);
        GR_MUST_SUCCEED(gr_set_fmpq(y, q, K));
        status = gr_tower_lazy_hurwitz_zeta(z, x, y, K);
        GR_MUST_SUCCEED(gr_mul(y, pi, pi, K));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 2, K));
        CHECK(status == GR_SUCCESS && gr_equal(z, y, K) == T_TRUE, "zeta(2, 1/2) = pi^2/2");

        /* zeta(3, 1/3): transcendental (conjecturally), checked numerically */
        GR_MUST_SUCCEED(gr_set_ui(x, 3, K));
        fmpq_set_si(q, 1, 3);
        GR_MUST_SUCCEED(gr_set_fmpq(y, q, K));
        status = gr_tower_lazy_hurwitz_zeta(z, x, y, K);
        if (status == GR_SUCCESS)
        {
            acb_set_ui(s, 3);
            acb_set_fmpq(a, q, PREC);
            acb_hurwitz_zeta(v, s, a, PREC);
            CHECK(_overlaps(z, v, K), "zeta(3, 1/3) numerically");
        }

        GR_TMP_CLEAR4(x, y, z, pi, K);
        fmpq_clear(q);
        acb_clear(s);
        acb_clear(a);
        acb_clear(v);
    }

    /* logarithms of integers with prime factors beyond the trial
       division limit (factored completely when word-sized, which
       40009 * 50021 < 2^32 is on all machines): the relation holds
       without a search, and the generators are the logs of the primes */
    {
        gr_ptr x, y, z;
        fmpz_t n;
        slong lev;
        gr_tower_struct * T;

        GR_TMP_INIT3(x, y, z, K);
        fmpz_init(n);
        fmpz_set_ui(n, 40009);
        fmpz_mul_ui(n, n, 50021);
        GR_MUST_SUCCEED(gr_set_fmpz(x, n, K));
        GR_MUST_SUCCEED(gr_log(x, x, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 40009, K));
        GR_MUST_SUCCEED(gr_log(y, y, K));
        GR_MUST_SUCCEED(gr_set_ui(z, 50021, K));
        GR_MUST_SUCCEED(gr_log(z, z, K));
        GR_MUST_SUCCEED(gr_add(y, y, z, K));
        CHECK(gr_equal(x, y, K) == T_TRUE, "log(pq) = log p + log q for large primes");
        T = gr_tower_lazy_get_tower(&lev, x, K);
        CHECK(T->num_trans == 2, "two logarithm generators");
        fmpz_clear(n);
        GR_TMP_CLEAR3(x, y, z, K);
    }

    /* conversions between contexts: an element of a second lazy context
       and of a real view, in both directions, and between views of
       different levels */
    {
        gr_ctx_t K2, R2;
        gr_ptr x, y, u, v;
        int status;

        gr_ctx_init_tower_lazy(K2, QQ, GR_TOWER_MERGE_EXPRESS);
        gr_ctx_init_tower_lazy_view(R2, K2, GR_TOWER_LAZY_REAL);
        GR_TMP_INIT2(x, y, K);
        GR_TMP_INIT(u, K2);
        GR_TMP_INIT(v, R2);

        /* x = sqrt(2) + sqrt(3) + pi in K */
        GR_MUST_SUCCEED(gr_set_ui(x, 2, K));
        GR_MUST_SUCCEED(gr_sqrt(x, x, K));
        GR_MUST_SUCCEED(gr_set_ui(y, 3, K));
        GR_MUST_SUCCEED(gr_sqrt(y, y, K));
        GR_MUST_SUCCEED(gr_add(x, x, y, K));
        GR_MUST_SUCCEED(gr_pi(y, K));
        GR_MUST_SUCCEED(gr_add(x, x, y, K));

        status = gr_set_other(u, x, K, K2);
        CHECK(status == GR_SUCCESS, "set_other between lazy contexts");
        status = gr_set_other(v, u, K2, R2);
        CHECK(status == GR_SUCCESS, "set_other into a real view of another context");
        status = gr_set_other(y, v, R2, K);
        CHECK(status == GR_SUCCESS && gr_equal(x, y, K) == T_TRUE, "round trip through another context and its view");

        /* a nonreal element is rejected by the real view of the other context */
        GR_MUST_SUCCEED(gr_i(y, K));
        GR_MUST_SUCCEED(gr_add(y, y, x, K));
        CHECK(gr_set_other(v, y, K, R2) == GR_DOMAIN, "nonreal element into a real view");

        GR_TMP_CLEAR2(x, y, K);
        GR_TMP_CLEAR(u, K2);
        GR_TMP_CLEAR(v, R2);
        gr_ctx_clear(R2);
        gr_ctx_clear(K2);
    }

    /* the public element structure and its accessors: representations,
       and polynomials in the first generator whatever the representation
       (in a fresh context: in K, sqrt(2) and zeta_5 may have become
       elements of larger fields made by the tests above) */
    {
        gr_ctx_t L, K3;
        gr_ptr x, y, u, v;
        fmpq_poly_t P, Q;
        fmpz_poly_t M, N;
        const char * exprs[] = {"exp(2*pi*i/5)^2 + 3/2*exp(2*pi*i/5) - 1/3",
                                "(1 + exp(2*pi*i/101)^50)^3 / 7",
                                "1 / (2 + 2^(1/3))"};
        slong i;

        gr_ctx_init_tower_lazy(L, QQ, 0);
        CHECK(gr_ctx_sizeof_elem(L) == sizeof(gr_tower_lazy_elem_struct), "element size");
        fmpq_poly_init(P);
        fmpq_poly_init(Q);
        fmpz_poly_init(M);
        fmpz_poly_init(N);
        GR_TMP_INIT2(x, y, L);

        /* (1 + sqrt(2))^3 = 7 + 5 sqrt(2), modulus x^2 - 2 */
        GR_MUST_SUCCEED(gr_set_ui(x, 2, L));
        GR_MUST_SUCCEED(gr_sqrt(x, x, L));
        GR_MUST_SUCCEED(gr_add_ui(x, x, 1, L));
        GR_MUST_SUCCEED(gr_pow_ui(x, x, 3, L));
        CHECK(gr_tower_lazy_repr(x, L) == GR_TOWER_LAZY_REPR_DENSE, "dense representation");
        CHECK(((gr_tower_lazy_elem_struct *) x)->repr == GR_TOWER_LAZY_REPR_DENSE, "repr field");
        GR_MUST_SUCCEED(gr_tower_lazy_get_fmpq_poly(P, M, x, L));
        fmpq_poly_set_str(Q, "2  7 5");
        fmpz_poly_set_str(N, "3  -2 0 1");
        CHECK(fmpq_poly_equal(P, Q) && fmpz_poly_equal(M, N), "(1 + sqrt(2))^3");
        CHECK(fmpq_poly_equal(&((gr_tower_lazy_elem_struct *) x)->elem.dense.poly, Q), "poly field");

        /* a rational number: a constant with modulus x */
        GR_MUST_SUCCEED(gr_set_si(y, -3, L));
        GR_MUST_SUCCEED(gr_div_ui(y, y, 7, L));
        CHECK(gr_tower_lazy_repr(y, L) == GR_TOWER_LAZY_REPR_RATIONAL, "rational representation");
        GR_MUST_SUCCEED(gr_tower_lazy_get_fmpq_poly(P, M, y, L));
        fmpq_poly_set_str(Q, "1  -3/7");
        fmpz_poly_set_str(N, "2  0 1");
        CHECK(fmpq_poly_equal(P, Q) && fmpz_poly_equal(M, N), "rational polynomial");

        /* two generators, or a transcendental one: not a polynomial in the first */
        GR_MUST_SUCCEED(gr_set_ui(y, 3, L));
        GR_MUST_SUCCEED(gr_sqrt(y, y, L));
        GR_MUST_SUCCEED(gr_add(y, y, x, L));
        CHECK(gr_tower_lazy_get_fmpq_poly(P, NULL, y, L) == GR_DOMAIN, "two generators");
        CHECK(gr_tower_lazy_repr(y, L) == GR_TOWER_LAZY_REPR_FLAT, "flat representation");
        GR_MUST_SUCCEED(gr_pi(y, L));
        CHECK(gr_tower_lazy_get_fmpq_poly(P, NULL, y, L) == GR_DOMAIN, "transcendental generator");

        /* the same coordinates without dense forms (the polynomial is
           canonical), and through a flat round trip of get_data */
        gr_ctx_init_tower_lazy(K3, QQ, GR_TOWER_MERGE_EXPRESS);
        GR_MUST_SUCCEED(gr_tower_lazy_ctx_set_option(K3, GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT, 0));
        GR_TMP_INIT(u, K3);
        GR_TMP_INIT(v, K3);
        for (i = 0; i < 3; i++)
        {
            GR_MUST_SUCCEED(gr_set_str(x, exprs[i], L));
            GR_MUST_SUCCEED(gr_mul_ui(x, x, 1, L));
            GR_MUST_SUCCEED(gr_set_str(u, exprs[i], K3));
            GR_MUST_SUCCEED(gr_mul_ui(u, u, 1, K3));
            CHECK(gr_tower_lazy_repr(x, L) == GR_TOWER_LAZY_REPR_DENSE, "dense form made");
            CHECK(gr_tower_lazy_repr(u, K3) == GR_TOWER_LAZY_REPR_FLAT, "no dense form when disabled");
            GR_MUST_SUCCEED(gr_tower_lazy_get_fmpq_poly(P, M, x, L));
            GR_MUST_SUCCEED(gr_tower_lazy_get_fmpq_poly(Q, N, u, K3));
            CHECK(fmpq_poly_equal(P, Q) && fmpz_poly_equal(M, N), "dense and flat coordinates");
            CHECK(fmpq_poly_length(P) < fmpz_poly_length(M), "reduced polynomial");
            /* get_data makes x flat in place; the polynomial is unchanged */
            (void) gr_tower_lazy_get_data(x, L);
            CHECK(gr_tower_lazy_repr(x, L) == GR_TOWER_LAZY_REPR_FLAT, "get_data makes the element flat");
            GR_MUST_SUCCEED(gr_tower_lazy_get_fmpq_poly(Q, NULL, x, L));
            CHECK(fmpq_poly_equal(P, Q), "coordinates after get_data");
            GR_MUST_SUCCEED(gr_set(v, u, K3));
            CHECK(gr_equal(v, u, K3) == T_TRUE, "copy");
        }
        GR_TMP_CLEAR(u, K3);
        GR_TMP_CLEAR(v, K3);
        gr_ctx_clear(K3);

        GR_TMP_CLEAR2(x, y, L);
        gr_ctx_clear(L);
        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
        fmpz_poly_clear(M);
        fmpz_poly_clear(N);
    }

    gr_ctx_clear(R);
    gr_ctx_clear(K);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
