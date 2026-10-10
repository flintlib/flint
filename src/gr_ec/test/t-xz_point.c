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
    The x-only ladder against the Weierstrass one.

    B v^2 = u^3 + A u^2 + u becomes y^2 = x^3 + a x + b under
    x = (u + A/3)/B, with a = (3 - A^2)/(3B^2) and
    b = (2A^3 - 9A)/(27B^3). Taking B = 1 loses nothing and keeps the map
    to x = u + A/3, so the two ladders can be compared coordinate by
    coordinate.
*/
static void
check_against_weierstrass(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ptr A, a24, a, b, u, x, y, t, three;
        gr_ec_xz_point_t P, Q, Q2, D;
        gr_ec_point_t W, WK;
        fmpz_t k, m, km;
        ulong p = n_randprime(state, 10 + n_randint(state, 12), 1);
        slong tries;
        int ok, have = 0;

        if (p <= 3 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        GR_TMP_INIT5(A, a24, a, b, u, R);
        GR_TMP_INIT4(x, y, t, three, R);

        /* a general B, so the whole transformation is exercised and not
           just the B = 1 special case */
        ok = (gr_set_ui(three, 3, R) == GR_SUCCESS)
            && (gr_set_ui(A, n_randint(state, p), R) == GR_SUCCESS)
            && (gr_set_ui(b, 1 + n_randint(state, p - 1), R) == GR_SUCCESS);

        ok = ok && (gr_ec_montgomery_a24(a24, A, R) == GR_SUCCESS);

        if (!ok || gr_ec_ctx_init_from_montgomery(E, R, A, b) != GR_SUCCESS)
        {
            GR_TMP_CLEAR4(x, y, t, three, R);
            GR_TMP_CLEAR5(A, a24, a, b, u, R);
            gr_ctx_clear(R);
            continue;
        }

        gr_ec_xz_point_init(P, R);
        gr_ec_xz_point_init(Q, R);
        gr_ec_xz_point_init(Q2, R);
        gr_ec_xz_point_init(D, R);
        gr_ec_point_init(W, E);
        gr_ec_point_init(WK, E);
        fmpz_init(k); fmpz_init(m); fmpz_init(km);

        /* an x whose Weierstrass image is a real point */
        for (tries = 0; tries < 60 && !have; tries++)
        {
            if (gr_set_ui(u, n_randint(state, p), R) != GR_SUCCESS)
                continue;

            if (gr_ec_montgomery_x_to_weierstrass(x, u, A, b, R) != GR_SUCCESS)
                continue;

            /* and the map must invert */
            FLINT_TEST(gr_ec_weierstrass_x_to_montgomery(t, x, A, b, R) == GR_SUCCESS);
            FLINT_TEST(gr_equal(t, u, R) == T_TRUE);

            if (gr_ec_point_lift_x(W, x, E) == GR_SUCCESS)
                have = 1;
        }

        if (!have)
            goto next;

        FLINT_TEST(gr_ec_xz_point_set_x(P, u, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_is_zero(P, R) == T_FALSE);

        fmpz_set_ui(k, 1 + n_randint(state, 1000));
        fmpz_set_ui(m, 1 + n_randint(state, 1000));
        fmpz_mul(km, k, m);

        /* the two ladders agree on the x-coordinate */
        FLINT_TEST(gr_ec_xz_point_mul_fmpz(Q, P, k, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_mul_fmpz(WK, W, k, E) == GR_SUCCESS);

        if (gr_ec_xz_point_is_zero(Q, R) == T_TRUE)
            FLINT_TEST(gr_ec_point_is_inf(WK, E) == T_TRUE);
        else
        {
            FLINT_TEST(gr_ec_point_is_inf(WK, E) == T_FALSE);
            FLINT_TEST(gr_ec_xz_point_get_x(t, Q, R) == GR_SUCCESS);
            FLINT_TEST(gr_ec_montgomery_x_to_weierstrass(y, t, A, b, R) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_get_affine(x, t, WK, E) == GR_SUCCESS);
            FLINT_TEST(gr_equal(y, x, R) == T_TRUE);
        }

        /* [m][k] P and [mk] P are the same x */
        FLINT_TEST(gr_ec_xz_point_mul_fmpz(Q2, Q, m, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_mul_fmpz(D, P, km, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_equal(Q2, D, R) != T_FALSE);

        /* doubling one way and the other */
        FLINT_TEST(gr_ec_xz_point_dbl(Q2, P, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_mul_ui(D, P, 2, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_equal(Q2, D, R) == T_TRUE);

        /* the differential addition, with the difference supplied */
        FLINT_TEST(gr_ec_xz_point_mul_ui(Q2, P, 3, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_mul_ui(D, P, 2, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_dadd(D, D, P, P, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_equal(Q2, D, R) == T_TRUE);

        /* zero, and aliasing */
        fmpz_zero(km);
        FLINT_TEST(gr_ec_xz_point_mul_fmpz(Q2, P, km, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_is_zero(Q2, R) == T_TRUE);

        FLINT_TEST(gr_ec_xz_point_set(Q2, P, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_dbl(Q2, Q2, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_dbl(D, P, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_equal(Q2, D, R) == T_TRUE);

next:
        fmpz_clear(k); fmpz_clear(m); fmpz_clear(km);
        gr_ec_point_clear(WK, E);
        gr_ec_point_clear(W, E);
        gr_ec_xz_point_clear(D, R);
        gr_ec_xz_point_clear(Q2, R);
        gr_ec_xz_point_clear(Q, R);
        gr_ec_xz_point_clear(P, R);
        gr_ec_ctx_clear(E);
        GR_TMP_CLEAR4(x, y, t, three, R);
        GR_TMP_CLEAR5(A, a24, a, b, u, R);
        gr_ctx_clear(R);
    }
}

/*
    Nothing in the x-only formulas tests whether a quantity vanishes, so
    they are meaningful over Z/n for a composite n, where the Weierstrass
    group law is not. The ladder must run there and still compose.
*/
static void
check_composite_modulus(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ptr A, a24, u;
        gr_ec_xz_point_t P, Q, Q2, D;
        fmpz_t k, m, km;
        ulong n = 6 + n_randint(state, 100000);

        if (gr_ctx_init_nmod(R, n) != GR_SUCCESS)
            continue;

        GR_TMP_INIT3(A, a24, u, R);

        if (gr_set_ui(A, n_randint(state, n), R) != GR_SUCCESS
                || gr_ec_montgomery_a24(a24, A, R) != GR_SUCCESS)
        {
            /* 4 is not invertible modulo an even n */
            GR_TMP_CLEAR3(A, a24, u, R);
            gr_ctx_clear(R);
            continue;
        }

        gr_ec_xz_point_init(P, R);
        gr_ec_xz_point_init(Q, R);
        gr_ec_xz_point_init(Q2, R);
        gr_ec_xz_point_init(D, R);
        fmpz_init(k); fmpz_init(m); fmpz_init(km);

        FLINT_TEST(gr_set_ui(u, n_randint(state, n), R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_set_x(P, u, R) == GR_SUCCESS);

        fmpz_set_ui(k, 1 + n_randint(state, 500));
        fmpz_set_ui(m, 1 + n_randint(state, 500));
        fmpz_mul(km, k, m);

        /* it runs, and the answer composes, with no notion of a group */
        FLINT_TEST(gr_ec_xz_point_mul_fmpz(Q, P, k, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_mul_fmpz(Q2, Q, m, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_mul_fmpz(D, P, km, a24, R) == GR_SUCCESS);
        FLINT_TEST(gr_ec_xz_point_equal(Q2, D, R) != T_FALSE);

        fmpz_clear(k); fmpz_clear(m); fmpz_clear(km);
        gr_ec_xz_point_clear(D, R);
        gr_ec_xz_point_clear(Q2, R);
        gr_ec_xz_point_clear(Q, R);
        gr_ec_xz_point_clear(P, R);
        GR_TMP_CLEAR3(A, a24, u, R);
        gr_ctx_clear(R);
    }
}

TEST_FUNCTION_START(gr_ec_xz_point, state)
{
    check_against_weierstrass(state);
    check_composite_modulus(state);

    TEST_FUNCTION_END(state);
}
