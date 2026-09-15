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

TEST_GR_FUNCTION_START(gr_ec_aff_point, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_aff_point_t P, Q, S, A, B;
        gr_ptr x, y;
        int status;

        gr_ec_test_ring(R, state);

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        gr_ec_aff_point_init(P, E);
        gr_ec_aff_point_init(Q, E);
        gr_ec_aff_point_init(S, E);
        gr_ec_aff_point_init(A, E);
        gr_ec_aff_point_init(B, E);
        GR_TMP_INIT2(x, y, R);

        FLINT_TEST(gr_ec_aff_point_is_inf(P, E) == T_TRUE);
        FLINT_TEST(gr_ec_aff_point_is_on_curve(P, E) == T_TRUE);
        FLINT_TEST(gr_ec_aff_point_equal(P, P, E) == T_TRUE);

        /* the identity has no affine coordinates */
        FLINT_TEST(gr_ec_aff_point_get_affine(x, y, P, E) == GR_DOMAIN);

        status = gr_ec_aff_point_randtest(P, state, E);
        status |= gr_ec_aff_point_randtest(Q, state, E);

        if (status != GR_SUCCESS)
        {
            if (status & GR_UNABLE)
                count_unable++;
            else
                count_domain++;

            goto cleanup;
        }

        count_success++;

        FLINT_TEST(gr_ec_aff_point_is_on_curve(P, E) != T_FALSE);
        FLINT_TEST(gr_ec_aff_point_is_on_curve(Q, E) != T_FALSE);

        FLINT_TEST(gr_ec_aff_point_set(A, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_aff_point_equal(A, P, E) == T_TRUE);

        FLINT_TEST(gr_ec_aff_point_set(B, Q, E) == GR_SUCCESS);
        gr_ec_aff_point_swap(A, B, E);
        FLINT_TEST(gr_ec_aff_point_equal(A, Q, E) == T_TRUE);
        FLINT_TEST(gr_ec_aff_point_equal(B, P, E) == T_TRUE);

        /* coordinates round trip */
        if (gr_ec_aff_point_get_affine(x, y, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_set_affine(A, x, y, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_aff_point_equal(A, P, E) == T_TRUE);
        }

        /* a point off the curve is rejected */
        if (gr_ec_aff_point_get_affine(x, y, P, E) == GR_SUCCESS
                && gr_add_ui(y, y, 1, R) == GR_SUCCESS)
        {
            int s2 = gr_ec_aff_point_set_affine(A, x, y, E);
            FLINT_TEST(s2 == GR_DOMAIN || s2 == GR_UNABLE || s2 == GR_SUCCESS);

            if (s2 == GR_SUCCESS)
                FLINT_TEST(gr_ec_aff_point_is_on_curve(A, E) != T_FALSE);
        }

        /* P + O = P */
        FLINT_TEST(gr_ec_aff_point_zero(A, E) == GR_SUCCESS);

        if (gr_ec_aff_point_add(B, P, A, E) == GR_SUCCESS)
            FLINT_TEST(gr_ec_aff_point_equal(B, P, E) != T_FALSE);

        /* P - P = O, -(-P) = P */
        if (gr_ec_aff_point_neg(A, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_is_on_curve(A, E) != T_FALSE);

            if (gr_ec_aff_point_add(B, P, A, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_aff_point_is_inf(B, E) != T_FALSE);

            if (gr_ec_aff_point_neg(B, A, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_aff_point_equal(B, P, E) != T_FALSE);
        }

        /* P + P = 2P */
        if (gr_ec_aff_point_add(A, P, P, E) == GR_SUCCESS
                && gr_ec_aff_point_dbl(B, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_equal(A, B, E) != T_FALSE);
            FLINT_TEST(gr_ec_aff_point_is_on_curve(B, E) != T_FALSE);
        }

        /* commutativity */
        if (gr_ec_aff_point_add(A, P, Q, E) == GR_SUCCESS
                && gr_ec_aff_point_add(B, Q, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_equal(A, B, E) == T_TRUE);
            FLINT_TEST(gr_ec_aff_point_is_on_curve(A, E) != T_FALSE);
        }

        /* associativity */
        if (gr_ec_aff_point_add(S, P, Q, E) == GR_SUCCESS
                && gr_ec_aff_point_add(A, S, Q, E) == GR_SUCCESS
                && gr_ec_aff_point_dbl(S, Q, E) == GR_SUCCESS
                && gr_ec_aff_point_add(B, P, S, E) == GR_SUCCESS)
            FLINT_TEST(gr_ec_aff_point_equal(A, B, E) == T_TRUE);

        /* P - Q = P + (-Q) */
        if (gr_ec_aff_point_sub(A, P, Q, E) == GR_SUCCESS
                && gr_ec_aff_point_neg(S, Q, E) == GR_SUCCESS
                && gr_ec_aff_point_add(B, P, S, E) == GR_SUCCESS)
            FLINT_TEST(gr_ec_aff_point_equal(A, B, E) == T_TRUE);

        /* aliasing */
        if (gr_ec_aff_point_add(A, P, Q, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_set(B, P, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_aff_point_add(B, B, Q, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_aff_point_equal(A, B, E) == T_TRUE);

            FLINT_TEST(gr_ec_aff_point_set(B, Q, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_aff_point_add(B, P, B, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_aff_point_equal(A, B, E) == T_TRUE);
        }

        if (gr_ec_aff_point_dbl(A, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_set(B, P, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_aff_point_dbl(B, B, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_aff_point_equal(A, B, E) == T_TRUE);
        }

        /* the group law needs a field */
        if (gr_ec_ctx_is_over_field(E) == T_FALSE)
            FLINT_TEST(gr_ec_aff_point_dbl(A, P, E) == GR_DOMAIN);

        /* on a short curve the general formulas must give the same
           result as the specialized ones */
        if (gr_ec_ctx_model(E) == GR_EC_SHORT_WEIERSTRASS
                && gr_ec_aff_point_is_inf(P, E) == T_FALSE
                && gr_ec_aff_point_is_inf(Q, E) == T_FALSE)
        {
            if (_gr_ec_aff_point_add_short_weierstrass(A, P, Q, E) == GR_SUCCESS
                    && _gr_ec_aff_point_add_long_weierstrass(B, P, Q, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_aff_point_equal(A, B, E) != T_FALSE);

            if (_gr_ec_aff_point_dbl_short_weierstrass(A, P, E) == GR_SUCCESS
                    && _gr_ec_aff_point_dbl_long_weierstrass(B, P, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_aff_point_equal(A, B, E) != T_FALSE);
        }

cleanup:
        GR_TMP_CLEAR2(x, y, R);
        gr_ec_aff_point_clear(P, E);
        gr_ec_aff_point_clear(Q, E);
        gr_ec_aff_point_clear(S, E);
        gr_ec_aff_point_clear(A, E);
        gr_ec_aff_point_clear(B, E);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
