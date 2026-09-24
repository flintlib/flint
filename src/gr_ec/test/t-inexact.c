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

/*
    Over a ring with inexact representation almost every predicate is
    allowed to come back T_UNKNOWN and almost every operation GR_UNABLE.
    What is not allowed is a definite wrong answer, which is what this
    test looks for.
*/
TEST_GR_FUNCTION_START(gr_ec_inexact, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, Q, A, B;
        gr_ec_jac_point_t Pj, Aj;
        gr_ec_aff_point_t Pa;
        int status;

        gr_ec_test_ring_inexact(R, state);

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        FLINT_TEST(gr_ec_ctx_is_smooth(E) != T_FALSE);

        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);
        gr_ec_point_init(A, E);
        gr_ec_point_init(B, E);
        gr_ec_jac_point_init(Pj, E);
        gr_ec_jac_point_init(Aj, E);
        gr_ec_aff_point_init(Pa, E);

        /* exact special values stay decidable */
        FLINT_TEST(gr_ec_point_is_inf(P, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_is_on_curve(P, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_equal(P, P, E) == T_TRUE);

        status = gr_ec_point_randtest(P, state, E);
        status |= gr_ec_point_randtest(Q, state, E);

        if (status != GR_SUCCESS)
        {
            if (status & GR_UNABLE)
                count_unable++;
            else
                count_domain++;

            goto cleanup;
        }

        count_success++;

        /* nothing below may ever come back as a definite negative */
        FLINT_TEST(gr_ec_point_is_on_curve(P, E) != T_FALSE);

        if (gr_ec_point_neg(A, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_point_is_on_curve(A, E) != T_FALSE);

            if (gr_ec_point_add(B, P, A, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_point_is_inf(B, E) != T_FALSE);
        }

        if (gr_ec_point_add(A, P, P, E) == GR_SUCCESS
                && gr_ec_point_dbl(B, P, E) == GR_SUCCESS)
            FLINT_TEST(gr_ec_point_equal(A, B, E) != T_FALSE);

        if (gr_ec_point_add(A, P, Q, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_point_is_on_curve(A, E) != T_FALSE);

            if (gr_ec_point_add(B, Q, P, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_point_equal(A, B, E) != T_FALSE);
        }

        /* conversions must not silently produce a valid but wrong point */
        if (gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_is_on_curve(Pj, E) != T_FALSE);

            if (gr_ec_point_set_jac_point(A, Pj, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_point_equal(A, P, E) != T_FALSE);

            if (gr_ec_jac_point_dbl(Aj, Pj, E) == GR_SUCCESS
                    && gr_ec_point_dbl(A, P, E) == GR_SUCCESS
                    && gr_ec_point_set_jac_point(B, Aj, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_point_equal(A, B, E) != T_FALSE);
        }
        else
        {
            FLINT_TEST(gr_ec_jac_point_is_inf(Pj, E) == T_UNKNOWN);
        }

        if (gr_ec_aff_point_set_point(Pa, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_is_on_curve(Pa, E) != T_FALSE);

            if (gr_ec_point_set_aff_point(A, Pa, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_point_equal(A, P, E) != T_FALSE);
        }

        /* scalar multiplication */
        if (gr_ec_point_mul_ui(A, P, 5, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_point_is_on_curve(A, E) != T_FALSE);

            {
                fmpz_t n;
                fmpz_init_set_ui(n, 5);

                if (_gr_ec_point_mul_fmpz_binary(B, P, n, E) == GR_SUCCESS)
                    FLINT_TEST(gr_ec_point_equal(A, B, E) != T_FALSE);

                fmpz_clear(n);
            }
        }

cleanup:
        gr_ec_aff_point_clear(Pa, E);
        gr_ec_jac_point_clear(Pj, E);
        gr_ec_jac_point_clear(Aj, E);
        gr_ec_point_clear(P, E);
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(A, E);
        gr_ec_point_clear(B, E);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
