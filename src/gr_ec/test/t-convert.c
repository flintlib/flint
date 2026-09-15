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

TEST_GR_FUNCTION_START(gr_ec_convert, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, Q, S, A;
        gr_ec_aff_point_t Pa, Qa, Sa;
        gr_ec_jac_point_t Pj, Qj, Sj;
        int status;

        gr_ec_test_ring(R, state);

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);
        gr_ec_point_init(S, E);
        gr_ec_point_init(A, E);
        gr_ec_aff_point_init(Pa, E);
        gr_ec_aff_point_init(Qa, E);
        gr_ec_aff_point_init(Sa, E);
        gr_ec_jac_point_init(Pj, E);
        gr_ec_jac_point_init(Qj, E);
        gr_ec_jac_point_init(Sj, E);

        /* the identity converts to the identity in every representation */
        FLINT_TEST(gr_ec_point_zero(P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_aff_point_set_point(Pa, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_aff_point_is_inf(Pa, E) == T_TRUE);
        FLINT_TEST(gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_jac_point_is_inf(Pj, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_set_aff_point(A, Pa, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_is_inf(A, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_set_jac_point(A, Pj, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_is_inf(A, E) == T_TRUE);
        FLINT_TEST(gr_ec_jac_point_set_aff_point(Pj, Pa, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_jac_point_is_inf(Pj, E) == T_TRUE);
        FLINT_TEST(gr_ec_aff_point_set_jac_point(Pa, Pj, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_aff_point_is_inf(Pa, E) == T_TRUE);

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

        /* projective -> Jacobian -> projective */
        FLINT_TEST(gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_jac_point_is_on_curve(Pj, E) != T_FALSE);
        FLINT_TEST(gr_ec_point_set_jac_point(A, Pj, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(A, P, E) != T_FALSE);

        /* projective -> affine -> projective */
        if (gr_ec_aff_point_set_point(Pa, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_is_on_curve(Pa, E) != T_FALSE);
            FLINT_TEST(gr_ec_point_set_aff_point(A, Pa, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_equal(A, P, E) != T_FALSE);

            /* affine -> Jacobian -> affine */
            FLINT_TEST(gr_ec_jac_point_set_aff_point(Pj, Pa, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_is_on_curve(Pj, E) != T_FALSE);

            if (gr_ec_aff_point_set_jac_point(Qa, Pj, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_aff_point_equal(Qa, Pa, E) == T_TRUE);
        }

        /* Jacobian -> affine -> Jacobian */
        FLINT_TEST(gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS);

        if (gr_ec_aff_point_set_jac_point(Pa, Pj, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_set_aff_point(Qj, Pa, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_equal(Qj, Pj, E) != T_FALSE);
        }

        /* the conversions are group homomorphisms */
        if (gr_ec_point_add(S, P, Q, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_set_point(Qj, Q, E) == GR_SUCCESS);

            if (gr_ec_jac_point_add(Sj, Pj, Qj, E) == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_set_jac_point(A, Sj, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(A, S, E) != T_FALSE);
            }

            if (gr_ec_aff_point_set_point(Pa, P, E) == GR_SUCCESS
                    && gr_ec_aff_point_set_point(Qa, Q, E) == GR_SUCCESS
                    && gr_ec_aff_point_add(Sa, Pa, Qa, E) == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_set_aff_point(A, Sa, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(A, S, E) != T_FALSE);

                /* mixed addition too */
                if (gr_ec_jac_point_add_aff_point(Sj, Pj, Qa, E) == GR_SUCCESS)
                {
                    FLINT_TEST(gr_ec_point_set_jac_point(A, Sj, E) == GR_SUCCESS);
                    FLINT_TEST(gr_ec_point_equal(A, S, E) != T_FALSE);
                }
            }
        }

        /* doubling commutes with the conversions */
        if (gr_ec_point_dbl(S, P, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS);

            if (gr_ec_jac_point_dbl(Sj, Pj, E) == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_set_jac_point(A, Sj, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(A, S, E) != T_FALSE);
            }

            if (gr_ec_aff_point_set_point(Pa, P, E) == GR_SUCCESS
                    && gr_ec_aff_point_dbl(Sa, Pa, E) == GR_SUCCESS)
            {
                FLINT_TEST(gr_ec_point_set_aff_point(A, Sa, E) == GR_SUCCESS);
                FLINT_TEST(gr_ec_point_equal(A, S, E) != T_FALSE);
            }
        }

cleanup:
        gr_ec_jac_point_clear(Pj, E);
        gr_ec_jac_point_clear(Qj, E);
        gr_ec_jac_point_clear(Sj, E);
        gr_ec_aff_point_clear(Pa, E);
        gr_ec_aff_point_clear(Qa, E);
        gr_ec_aff_point_clear(Sa, E);
        gr_ec_point_clear(P, E);
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(S, E);
        gr_ec_point_clear(A, E);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
