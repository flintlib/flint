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
#include "t-helpers.h"

/*
    A fresh context knows nothing, asking for the order counts and
    remembers, and what comes back afterwards is the same number.
*/
static void
check_cache(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t N, M;
        ulong p = n_randprime(state, 5 + n_randint(state, 10), 1);

        if (p <= 3 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(N);
        fmpz_init(M);

        FLINT_TEST(gr_ec_ctx_order_kind(E) == GR_EC_ORDER_UNKNOWN);
        FLINT_TEST(gr_ec_ctx_get_cached_order(M, E) == GR_EC_ORDER_UNKNOWN);

        if (gr_ec_ctx_order(N, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_ctx_order_kind(E) == GR_EC_ORDER_EXACT);
            FLINT_TEST(gr_ec_ctx_get_cached_order(M, E) == GR_EC_ORDER_EXACT);
            FLINT_TEST(fmpz_equal(M, N));

            /* and it agrees with counting the hard way */
            if (gr_ec_ctx_cardinality_naive(M, E) == GR_SUCCESS)
                FLINT_TEST(fmpz_equal(M, N));

            /* forgetting it puts the context back where it started */
            gr_ec_ctx_clear_order(E);
            FLINT_TEST(gr_ec_ctx_order_kind(E) == GR_EC_ORDER_UNKNOWN);
        }

        fmpz_clear(N);
        fmpz_clear(M);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/* A hint is accepted when it is right and refused when it is not. */
static void
check_hint(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t N, M;
        ulong p = n_randprime(state, 5 + n_randint(state, 10), 1);

        if (p <= 3 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(N);
        fmpz_init(M);

        if (gr_ec_ctx_cardinality_naive(N, E) != GR_SUCCESS)
            goto next;

        /* the truth, on a context that has not been told anything */
        FLINT_TEST(gr_ec_ctx_set_order(E, N) == GR_SUCCESS);
        FLINT_TEST(gr_ec_ctx_order_kind(E) == GR_EC_ORDER_EXACT);
        FLINT_TEST(gr_ec_ctx_get_cached_order(M, E) == GR_EC_ORDER_EXACT);
        FLINT_TEST(fmpz_equal(M, N));

        /* anything else is refused, now that the right answer is known */
        fmpz_add_ui(M, N, 1);
        FLINT_TEST(gr_ec_ctx_set_order(E, M) == GR_DOMAIN);
        fmpz_set_si(M, -1);
        FLINT_TEST(gr_ec_ctx_set_order(E, M) == GR_DOMAIN);
        fmpz_zero(M);
        FLINT_TEST(gr_ec_ctx_set_order(E, M) == GR_DOMAIN);

        /* the cache still holds the right value */
        FLINT_TEST(gr_ec_ctx_get_cached_order(M, E) == GR_EC_ORDER_EXACT);
        FLINT_TEST(fmpz_equal(M, N));

        /* multiples of the order annihilate; N + 1 generally does not */
        FLINT_TEST(gr_ec_ctx_annihilates(E, N) == T_TRUE);
        fmpz_mul_ui(M, N, 3);
        FLINT_TEST(gr_ec_ctx_annihilates(E, M) == T_TRUE);

next:
        fmpz_clear(N);
        fmpz_clear(M);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    A wrong hint on a context that has not counted anything has to be caught
    by the random point probe rather than by a comparison, so drive that
    path on purpose with a field too big for the exact check.
*/
static void
check_hint_probe(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t N, M;
        ulong p = n_randprime(state, 20 + n_randint(state, 6), 1);

        if (gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(N);
        fmpz_init(M);

        if (gr_ec_ctx_cardinality_bsgs(N, E) != GR_SUCCESS)
            goto next;

        /* a value inside the Hasse interval but not the order: only the
           probe can tell, and it has to */
        fmpz_add_ui(M, N, 1);

        if (gr_ec_ctx_annihilates(E, M) == T_FALSE)
            FLINT_TEST(gr_ec_ctx_set_order(E, M) == GR_DOMAIN);

        FLINT_TEST(gr_ec_ctx_order_kind(E) == GR_EC_ORDER_UNKNOWN);

        /* the truth passes the same probe */
        FLINT_TEST(gr_ec_ctx_set_order(E, N) == GR_SUCCESS);

        /* and something wildly outside the Hasse interval never gets that far */
        fmpz_mul_ui(M, N, 100);
        fmpz_add_ui(M, M, 1);
        gr_ec_ctx_clear_order(E);
        FLINT_TEST(gr_ec_ctx_set_order(E, M) == GR_DOMAIN);

next:
        fmpz_clear(N);
        fmpz_clear(M);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    Scalars taken from Z/n: the ring acts exactly when n kills the group,
    and the answer is what multiplying by a lift of the scalar gives.
*/
static void
check_scalar_ring(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R, S;
        gr_ec_ctx_t E;
        gr_ec_point_t P, A, B;
        fmpz_t N, k;
        gr_ptr s;
        ulong p = n_randprime(state, 5 + n_randint(state, 10), 1);

        if (p <= 3 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(N);
        fmpz_init(k);
        gr_ec_point_init(P, E);
        gr_ec_point_init(A, E);
        gr_ec_point_init(B, E);

        if (gr_ec_ctx_order(N, E) != GR_SUCCESS
                || !fmpz_abs_fits_ui(N) || fmpz_cmp_ui(N, 2) < 0
                || gr_ec_point_randtest(P, state, E) != GR_SUCCESS)
            goto next;

        if (gr_ctx_init_nmod(S, fmpz_get_ui(N)) == GR_SUCCESS)
        {
            GR_TMP_INIT(s, S);

            if (gr_randtest(s, state, S) == GR_SUCCESS
                    && gr_get_fmpz(k, s, S) == GR_SUCCESS)
            {
                FLINT_TEST(gr_mul_other(A, P, s, S, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_mul_fmpz(B, P, k, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);

                /* the same thing written the other way round */
                FLINT_TEST(gr_other_mul(A, s, S, P, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);
            }

            GR_TMP_CLEAR(s, S);
            gr_ctx_clear(S);
        }

        /* a modulus that leaves some point alive is refused outright */
        {
            fmpz_t bad;
            fmpz_init(bad);
            fmpz_add_ui(bad, N, 1);

            if (fmpz_abs_fits_ui(bad) && fmpz_cmp_ui(bad, 2) >= 0
                    && gr_ec_ctx_annihilates(E, bad) == T_FALSE
                    && gr_ctx_init_nmod(S, fmpz_get_ui(bad)) == GR_SUCCESS)
            {
                GR_TMP_INIT(s, S);

                if (gr_randtest(s, state, S) == GR_SUCCESS)
                    FLINT_TEST(gr_mul_other(A, P, s, S, E) == GR_DOMAIN);

                GR_TMP_CLEAR(s, S);
                gr_ctx_clear(S);
            }

            fmpz_clear(bad);
        }

        /* and a ring that is not a Z/n at all never acts */
        {
            gr_ctx_t Q;
            gr_ptr z;

            gr_ctx_init_fmpq(Q);
            GR_TMP_INIT(z, Q);
            FLINT_TEST(gr_one(z, Q) == GR_SUCCESS);
            FLINT_TEST(gr_mul_other(A, P, z, Q, E) == GR_DOMAIN);
            GR_TMP_CLEAR(z, Q);
            gr_ctx_clear(Q);
        }

next:
        gr_ec_point_clear(B, E);
        gr_ec_point_clear(A, E);
        gr_ec_point_clear(P, E);
        fmpz_clear(N);
        fmpz_clear(k);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    Reducing the scalar against a cached annihilator must not change any
    answer, however large the scalar is.
*/
static void
check_reduction(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E, F;
        gr_ec_point_t P, A, B;
        fmpz_t N, big;
        ulong p = n_randprime(state, 6 + n_randint(state, 10), 1);

        if (p <= 3 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(N);
        fmpz_init(big);
        gr_ec_point_init(P, E);
        gr_ec_point_init(A, E);
        gr_ec_point_init(B, E);

        if (gr_ec_point_randtest(P, state, E) != GR_SUCCESS)
            goto next;

        fmpz_randbits(big, state, 1 + n_randint(state, 200));

        /* with nothing cached */
        if (gr_ec_point_mul_fmpz(A, P, big, E) != GR_SUCCESS)
            goto next;

        /* and with the order cached, which turns on the reduction */
        if (gr_ec_ctx_order(N, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_point_mul_fmpz(B, P, big, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);
        }

        /* an annihilator that is not the order reduces just as validly */
        if (gr_ec_ctx_order(N, E) == GR_SUCCESS
                && gr_ec_ctx_init_randtest(F, state, R) == GR_SUCCESS)
        {
            gr_ec_ctx_clear(F);

            fmpz_mul_ui(N, N, 6);
            gr_ec_ctx_clear_order(E);
            FLINT_TEST(gr_ec_ctx_set_annihilator(E, N) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_mul_fmpz(B, P, big, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);
        }

next:
        gr_ec_point_clear(B, E);
        gr_ec_point_clear(A, E);
        gr_ec_point_clear(P, E);
        fmpz_clear(N);
        fmpz_clear(big);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

TEST_FUNCTION_START(gr_ec_order, state)
{
    check_cache(state);
    check_hint(state);
    check_hint_probe(state);
    check_scalar_ring(state);
    check_reduction(state);

    TEST_FUNCTION_END(state);
}
