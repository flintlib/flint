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

TEST_GR_FUNCTION_START(gr_ec_jac_point, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_jac_point_t P, Q, S, A, B;
        gr_ec_aff_point_t Qa;
        gr_ptr x, y;
        int status;

        gr_ec_test_ring(R, state);

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        gr_ec_jac_point_init(P, E);
        gr_ec_jac_point_init(Q, E);
        gr_ec_jac_point_init(S, E);
        gr_ec_jac_point_init(A, E);
        gr_ec_jac_point_init(B, E);
        gr_ec_aff_point_init(Qa, E);
        GR_TMP_INIT2(x, y, R);

        FLINT_TEST(gr_ec_jac_point_is_inf(P, E) == T_TRUE);
        FLINT_TEST(gr_ec_jac_point_is_on_curve(P, E) == T_TRUE);
        FLINT_TEST(gr_ec_jac_point_equal(P, P, E) == T_TRUE);
        FLINT_TEST(gr_ec_jac_point_get_affine(x, y, P, E) == GR_DOMAIN);

        status = gr_ec_jac_point_randtest(P, state, E);
        status |= gr_ec_jac_point_randtest(Q, state, E);

        if (status != GR_SUCCESS)
        {
            if (status & GR_UNABLE)
                count_unable++;
            else
                count_domain++;

            goto cleanup;
        }

        count_success++;

        FLINT_TEST(gr_ec_jac_point_is_on_curve(P, E) != T_FALSE);
        FLINT_TEST(gr_ec_jac_point_is_on_curve(Q, E) != T_FALSE);

        FLINT_TEST(gr_ec_jac_point_set(A, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_jac_point_equal(A, P, E) == T_TRUE);

        FLINT_TEST(gr_ec_jac_point_set(B, Q, E) == GR_SUCCESS);
        gr_ec_jac_point_swap(A, B, E);
        FLINT_TEST(gr_ec_jac_point_equal(A, Q, E) == T_TRUE);
        FLINT_TEST(gr_ec_jac_point_equal(B, P, E) == T_TRUE);

        /* normalization preserves the point and sets Z = 1 */
        if (gr_ec_jac_point_normalize(A, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_equal(A, P, E) != T_FALSE);
            FLINT_TEST(gr_ec_jac_point_is_on_curve(A, E) != T_FALSE);

            if (gr_ec_jac_point_is_inf(A, E) == T_FALSE)
                FLINT_TEST(gr_is_one(gr_ec_jac_point_z_srcptr(A, E), R) != T_FALSE);
        }

        /* coordinates round trip */
        if (gr_ec_jac_point_get_affine(x, y, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_set_affine(A, x, y, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_equal(A, P, E) == T_TRUE);
        }

        /* raw Jacobian coordinates round trip */
        if (gr_ec_jac_point_is_inf(P, E) == T_FALSE)
        {
            status = gr_ec_jac_point_set_jacobian(A,
                        gr_ec_jac_point_x_srcptr(P, E),
                        gr_ec_jac_point_y_srcptr(P, E),
                        gr_ec_jac_point_z_srcptr(P, E), E);

            if (status == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_equal(A, P, E) == T_TRUE);
        }

        /* P + O = P */
        FLINT_TEST(gr_ec_jac_point_zero(A, E) == GR_SUCCESS);

        if (gr_ec_jac_point_add(B, P, A, E) == GR_SUCCESS)
            FLINT_TEST(gr_ec_jac_point_equal(B, P, E) != T_FALSE);

        /* P - P = O, -(-P) = P */
        if (gr_ec_jac_point_neg(A, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_is_on_curve(A, E) != T_FALSE);

            if (gr_ec_jac_point_add(B, P, A, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_is_inf(B, E) != T_FALSE);

            if (gr_ec_jac_point_neg(B, A, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_equal(B, P, E) != T_FALSE);
        }

        /* P + P = 2P */
        if (gr_ec_jac_point_add(A, P, P, E) == GR_SUCCESS
                && gr_ec_jac_point_dbl(B, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);
            FLINT_TEST(gr_ec_jac_point_is_on_curve(B, E) != T_FALSE);
        }

        /* commutativity and associativity */
        if (gr_ec_jac_point_add(A, P, Q, E) == GR_SUCCESS
                && gr_ec_jac_point_add(B, Q, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);
            FLINT_TEST(gr_ec_jac_point_is_on_curve(A, E) != T_FALSE);
        }

        if (gr_ec_jac_point_add(S, P, Q, E) == GR_SUCCESS
                && gr_ec_jac_point_add(A, S, Q, E) == GR_SUCCESS
                && gr_ec_jac_point_dbl(S, Q, E) == GR_SUCCESS
                && gr_ec_jac_point_add(B, P, S, E) == GR_SUCCESS)
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);

        /* P - Q = P + (-Q) */
        if (gr_ec_jac_point_sub(A, P, Q, E) == GR_SUCCESS
                && gr_ec_jac_point_neg(S, Q, E) == GR_SUCCESS
                && gr_ec_jac_point_add(B, P, S, E) == GR_SUCCESS)
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);

        /* mixed addition agrees with the general one */
        if (gr_ec_aff_point_set_jac_point(Qa, Q, E) == GR_SUCCESS)
        {
            if (gr_ec_jac_point_add(A, P, Q, E) == GR_SUCCESS
                    && gr_ec_jac_point_add_aff_point(B, P, Qa, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);

            if (gr_ec_jac_point_sub(A, P, Q, E) == GR_SUCCESS
                    && gr_ec_jac_point_sub_aff_point(B, P, Qa, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);
        }

        /* aliasing */
        if (gr_ec_jac_point_add(A, P, Q, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_set(B, P, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_add(B, B, Q, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);

            FLINT_TEST(gr_ec_jac_point_set(B, Q, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_add(B, P, B, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);
        }

        if (gr_ec_jac_point_dbl(A, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_set(B, P, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_dbl(B, B, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);
        }

        /* on a short curve the general formulas must give the same
           result as the specialized ones */
        if (gr_ec_ctx_model(E) == GR_EC_SHORT_WEIERSTRASS
                && gr_ec_jac_point_is_inf(P, E) == T_FALSE
                && gr_ec_jac_point_is_inf(Q, E) == T_FALSE)
        {
            if (_gr_ec_jac_point_add_short_weierstrass(A, P, Q, E) == GR_SUCCESS
                    && _gr_ec_jac_point_add_long_weierstrass(B, P, Q, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);

            if (_gr_ec_jac_point_dbl_short_weierstrass(A, P, E) == GR_SUCCESS
                    && _gr_ec_jac_point_dbl_long_weierstrass(B, P, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);

            if (gr_ec_aff_point_set_jac_point(Qa, Q, E) == GR_SUCCESS
                    && _gr_ec_jac_point_add_aff_point_short_weierstrass(A, P, Qa, E) == GR_SUCCESS
                    && _gr_ec_jac_point_add_aff_point_long_weierstrass(B, P, Qa, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE);
        }

cleanup:
        GR_TMP_CLEAR2(x, y, R);
        gr_ec_aff_point_clear(Qa, E);
        gr_ec_jac_point_clear(P, E);
        gr_ec_jac_point_clear(Q, E);
        gr_ec_jac_point_clear(S, E);
        gr_ec_jac_point_clear(A, E);
        gr_ec_jac_point_clear(B, E);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
