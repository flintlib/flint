/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "t-helpers.h"

/* curves over Q with known invariants */
static void
check_known_curves(void)
{
    gr_ctx_t R;
    gr_ec_ctx_t E;
    gr_ptr t, u;

    gr_ctx_init_fmpq(R);
    GR_TMP_INIT2(t, u, R);

    /* y^2 = x^3 - x : disc = 64, j = 1728 */
    FLINT_TEST(gr_ec_ctx_init_short_weierstrass_si(E, R, -1, 0) == GR_SUCCESS);
    FLINT_TEST(gr_ec_ctx_model(E) == GR_EC_SHORT_WEIERSTRASS);
    FLINT_TEST(gr_ec_ctx_discriminant(t, E) == GR_SUCCESS);
    FLINT_TEST(gr_set_si(u, 64, R) == GR_SUCCESS);
    FLINT_TEST(gr_equal(t, u, R) == T_TRUE);
    FLINT_TEST(gr_ec_ctx_j_invariant(t, E) == GR_SUCCESS);
    FLINT_TEST(gr_set_si(u, 1728, R) == GR_SUCCESS);
    FLINT_TEST(gr_equal(t, u, R) == T_TRUE);
    gr_ec_ctx_clear(E);

    /* y^2 = x^3 + 1 : disc = -432, j = 0 */
    FLINT_TEST(gr_ec_ctx_init_short_weierstrass_si(E, R, 0, 1) == GR_SUCCESS);
    FLINT_TEST(gr_ec_ctx_discriminant(t, E) == GR_SUCCESS);
    FLINT_TEST(gr_set_si(u, -432, R) == GR_SUCCESS);
    FLINT_TEST(gr_equal(t, u, R) == T_TRUE);
    FLINT_TEST(gr_ec_ctx_j_invariant(t, E) == GR_SUCCESS);
    FLINT_TEST(gr_is_zero(t, R) == T_TRUE);
    gr_ec_ctx_clear(E);

    /* y^2 + y = x^3 - x^2 - 10x - 20 (the curve 11a1): disc = -161051, j = -122023936/161051 */
    FLINT_TEST(gr_ec_ctx_init_si(E, R, 0, -1, 1, -10, -20) == GR_SUCCESS);
    FLINT_TEST(gr_ec_ctx_model(E) == GR_EC_LONG_WEIERSTRASS);
    FLINT_TEST(gr_ec_ctx_discriminant(t, E) == GR_SUCCESS);
    FLINT_TEST(gr_set_si(u, -161051, R) == GR_SUCCESS);
    FLINT_TEST(gr_equal(t, u, R) == T_TRUE);
    gr_ec_ctx_clear(E);

    /* singular curves are rejected */
    FLINT_TEST(gr_ec_ctx_init_short_weierstrass_si(E, R, 0, 0) == GR_DOMAIN);
    FLINT_TEST(gr_ec_ctx_init_si(E, R, 0, 0, 0, 0, 0) == GR_DOMAIN);

    GR_TMP_CLEAR2(t, u, R);
    gr_ctx_clear(R);
}

TEST_FUNCTION_START(gr_ec_ctx, state)
{
    slong iter;

    check_known_curves();

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ptr b2, b4, b6, b8, c4, c6, t, u;
        int status;

        gr_ec_test_ring(R, state);

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        /* a context is never created with a provably zero discriminant */
        FLINT_TEST(gr_ec_ctx_is_smooth(E) != T_FALSE);

        /* the model must match the stored a-invariants */
        if (gr_ec_ctx_model(E) == GR_EC_SHORT_WEIERSTRASS)
        {
            FLINT_TEST(gr_is_zero(GR_EC_A1(E), R) == T_TRUE);
            FLINT_TEST(gr_is_zero(GR_EC_A2(E), R) == T_TRUE);
            FLINT_TEST(gr_is_zero(GR_EC_A3(E), R) == T_TRUE);
        }

        GR_TMP_INIT5(b2, b4, b6, b8, c4, R);
        GR_TMP_INIT3(c6, t, u, R);

        status = gr_ec_ctx_b_invariants(b2, b4, b6, b8, E);
        status |= gr_ec_ctx_c_invariants(c4, c6, E);

        /* 4 b8 = b2 b6 - b4^2 */
        if (status == GR_SUCCESS)
        {
            status |= gr_mul_ui(t, b8, 4, R);
            status |= gr_mul(u, b2, b6, R);
            status |= gr_submul(u, b4, b4, R);
            status |= gr_sub(t, t, u, R);

            if (status == GR_SUCCESS)
                FLINT_TEST(gr_is_zero(t, R) != T_FALSE);
        }

        /* c4^3 - c6^2 = 1728 disc */
        status = gr_sqr(t, c4, R);
        status |= gr_mul(t, t, c4, R);
        status |= gr_submul(t, c6, c6, R);
        status |= gr_ec_ctx_discriminant(u, E);
        status |= gr_mul_ui(u, u, 1728, R);
        status |= gr_sub(t, t, u, R);

        if (status == GR_SUCCESS)
            FLINT_TEST(gr_is_zero(t, R) != T_FALSE);

        /* j disc = c4^3 whenever j exists */
        if (gr_ec_ctx_j_invariant(t, E) == GR_SUCCESS)
        {
            status = gr_ec_ctx_discriminant(u, E);
            status |= gr_mul(t, t, u, R);
            status |= gr_sqr(u, c4, R);
            status |= gr_mul(u, u, c4, R);
            status |= gr_sub(t, t, u, R);

            if (status == GR_SUCCESS)
                FLINT_TEST(gr_is_zero(t, R) != T_FALSE);
        }

        GR_TMP_CLEAR3(c6, t, u, R);
        GR_TMP_CLEAR5(b2, b4, b6, b8, c4, R);

        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}
