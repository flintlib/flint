/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpq.h"
#include "t-helpers.h"

#define MAX_POINTS 600

/*
    Every point of E(F_p), by walking x over the field and lifting. Only
    usable over a small prime field, which is the point: it gives an
    independent answer to "how many Q have n Q = P", against which the
    division functions can be held.
*/
static slong
all_points(gr_ec_point_struct * pts, slong maxpts, ulong p, gr_ec_ctx_t E)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(E);
    gr_ec_point_t Q, Qn;
    gr_ptr x;
    slong n = 0, i;

    gr_ec_point_init(Q, E);
    gr_ec_point_init(Qn, E);
    GR_TMP_INIT(x, R);

    gr_ec_point_init(pts + n, E);

    if (gr_ec_point_zero(pts + n, E) == GR_SUCCESS)
        n++;

    for (i = 0; i < (slong) p && n + 2 <= maxpts; i++)
    {
        if (gr_set_ui(x, (ulong) i, R) != GR_SUCCESS)
            continue;

        if (gr_ec_point_lift_x(Q, x, E) != GR_SUCCESS)
            continue;

        gr_ec_point_init(pts + n, E);

        if (gr_ec_point_set(pts + n, Q, E) == GR_SUCCESS)
            n++;

        if (gr_ec_point_neg(Qn, Q, E) == GR_SUCCESS
                && gr_ec_point_equal(Qn, Q, E) != T_TRUE)
        {
            gr_ec_point_init(pts + n, E);

            if (gr_ec_point_set(pts + n, Qn, E) == GR_SUCCESS)
                n++;
        }
    }

    GR_TMP_CLEAR(x, R);
    gr_ec_point_clear(Qn, E);
    gr_ec_point_clear(Q, E);

    return n;
}

static void
check_against_brute_force(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 400 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_struct pts[MAX_POINTS];
        gr_ec_point_t P, Q, T;
        fmpz_t n, N;
        slong npts = 0, i, sols = 0;
        ulong p = n_randprime(state, 4 + n_randint(state, 5), 1);
        ulong nn;
        int st;

        if (p <= 3 || p > 250 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        /* the root finding path is written for the short model */
        if (gr_ec_ctx_model(E) != GR_EC_SHORT_WEIERSTRASS)
        {
            gr_ec_ctx_clear(E);
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(n);
        fmpz_init(N);
        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);
        gr_ec_point_init(T, E);

        npts = all_points(pts, MAX_POINTS, p, E);

        /* the enumeration and the point count must be the same number */
        if (gr_ec_ctx_order(N, E) == GR_SUCCESS)
            FLINT_TEST(fmpz_equal_si(N, npts));

        if (gr_ec_point_randtest(P, state, E) != GR_SUCCESS)
            goto next;

        nn = 1 + n_randint(state, 12);
        fmpz_set_ui(n, nn);

        for (i = 0; i < npts; i++)
            if (gr_ec_point_mul_ui(T, pts + i, nn, E) == GR_SUCCESS
                    && gr_ec_point_equal(T, P, E) == T_TRUE)
                sols++;

        /* div wants exactly one */
        st = gr_ec_point_div_fmpz(Q, P, n, E);

        if (sols == 1)
        {
            FLINT_TEST(st == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_mul_ui(T, Q, nn, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_equal(T, P, E) == T_TRUE);
        }
        else
            FLINT_TEST(st == GR_DOMAIN || st == GR_UNABLE);

        /* div_nonunique wants at least one */
        st = gr_ec_point_div_fmpz_nonunique(Q, P, n, E);

        if (sols >= 1)
        {
            FLINT_TEST(st == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_mul_ui(T, Q, nn, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_equal(T, P, E) == T_TRUE);
        }
        else
            FLINT_TEST(st == GR_DOMAIN || st == GR_UNABLE);

        /* the ui and si spellings agree with the fmpz one */
        {
            gr_ec_point_t U;
            gr_ec_point_init(U, E);

            FLINT_TEST(gr_ec_point_div_ui(U, P, nn, E)
                    == gr_ec_point_div_fmpz(Q, P, n, E));

            if (gr_ec_point_div_fmpz(Q, P, n, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_point_equal(U, Q, E) == T_TRUE);

            /* dividing by -n is negating the answer */
            fmpz_neg(n, n);

            if (gr_ec_point_div_fmpz(U, P, n, E) == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_neg(T, U, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_div_si(U, P, (slong) nn, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(T, U, E) == T_TRUE);
            }

            fmpz_neg(n, n);
            gr_ec_point_clear(U, E);
        }

        /* (a/b) P is the Q with b Q = a P */
        {
            fmpq_t c;
            gr_ec_point_t U;

            fmpq_init(c);
            gr_ec_point_init(U, E);

            fmpq_set_si(c, 1 + (slong) n_randint(state, 9),
                    1 + n_randint(state, 9));

            if (gr_ec_point_mul_fmpq_nonunique(Q, P, c, E) == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_mul_fmpz(T, Q, fmpq_denref(c), E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_mul_fmpz(U, P, fmpq_numref(c), E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(T, U, E) == T_TRUE);
            }

            gr_ec_point_clear(U, E);
            fmpq_clear(c);
        }

        /* dividing by zero is never a single answer */
        fmpz_zero(n);
        FLINT_TEST(gr_ec_point_div_fmpz(Q, P, n, E) == GR_DOMAIN);

next:
        for (i = 0; i < npts; i++)
            gr_ec_point_clear(pts + i, E);

        gr_ec_point_clear(T, E);
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(P, E);
        fmpz_clear(n);
        fmpz_clear(N);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    The three representations must give the same answers, and the generic
    interface must be the module functions.
*/
static void
check_representations(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, Q;
        gr_ec_aff_point_t Pa, Qa;
        gr_ec_jac_point_t Pj, Qj;
        gr_ec_point_t Qc;
        fmpz_t n;
        ulong p = n_randprime(state, 5 + n_randint(state, 8), 1);
        ulong nn;
        int sp;

        if (p <= 3 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(n);
        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);
        gr_ec_point_init(Qc, E);
        gr_ec_aff_point_init(Pa, E);
        gr_ec_aff_point_init(Qa, E);
        gr_ec_jac_point_init(Pj, E);
        gr_ec_jac_point_init(Qj, E);

        nn = 1 + n_randint(state, 8);
        fmpz_set_ui(n, nn);

        if (gr_ec_point_randtest(P, state, E) != GR_SUCCESS)
            goto next;

        sp = gr_ec_point_div_fmpz(Q, P, n, E);

        if (gr_ec_aff_point_set_point(Pa, P, E) == GR_SUCCESS)
        {
            int sa = gr_ec_aff_point_div_fmpz(Qa, Pa, n, E);

            if (sp == GR_SUCCESS && sa == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_set_aff_point(Qc, Qa, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(Qc, Q, E) == T_TRUE);
            }
            else
                FLINT_TEST(sp == sa || sa == GR_UNABLE);
        }

        if (gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS)
        {
            int sj = gr_ec_jac_point_div_fmpz(Qj, Pj, n, E);

            if (sp == GR_SUCCESS && sj == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_set_jac_point(Qc, Qj, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(Qc, Q, E) == T_TRUE);
            }
            else
                FLINT_TEST(sp == sj || sj == GR_UNABLE);
        }

        /* the generic interface reaches the same code */
        FLINT_TEST(gr_div_fmpz(Qc, P, n, E) == sp);

        if (sp == GR_SUCCESS)
            FLINT_TEST(gr_ec_point_equal(Qc, Q, E) == T_TRUE);

next:
        gr_ec_jac_point_clear(Qj, E);
        gr_ec_jac_point_clear(Pj, E);
        gr_ec_aff_point_clear(Qa, E);
        gr_ec_aff_point_clear(Pa, E);
        gr_ec_point_clear(Qc, E);
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(P, E);
        fmpz_clear(n);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    Over a field big enough that the group order is coprime to small n, the
    cheap path is the only one taken, and it has to be right there too.
*/
static void
check_coprime_path(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, Q, T;
        fmpz_t n, N, g;
        ulong p = n_randprime(state, 20 + n_randint(state, 8), 1);
        ulong nn;

        if (gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(n);
        fmpz_init(N);
        fmpz_init(g);
        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);
        gr_ec_point_init(T, E);

        if (gr_ec_ctx_order(N, E) != GR_SUCCESS
                || gr_ec_point_randtest(P, state, E) != GR_SUCCESS)
            goto next;

        nn = 1 + n_randint(state, 30);
        fmpz_set_ui(n, nn);
        fmpz_gcd(g, n, N);

        if (!fmpz_is_one(g))
            goto next;

        /* coprime to the order, so there is exactly one answer */
        FLINT_TEST(gr_ec_point_div_fmpz(Q, P, n, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_mul_ui(T, Q, nn, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(T, P, E) == T_TRUE);

next:
        gr_ec_point_clear(T, E);
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(P, E);
        fmpz_clear(n);
        fmpz_clear(N);
        fmpz_clear(g);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

TEST_FUNCTION_START(gr_ec_div, state)
{
    check_against_brute_force(state);
    check_representations(state);
    check_coprime_path(state);

    TEST_FUNCTION_END(state);
}
