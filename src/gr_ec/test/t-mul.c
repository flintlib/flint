/*
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Opus 5 and Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "t-helpers.h"

TEST_GR_FUNCTION_START(gr_ec_mul, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_jac_point_t P, A, B, C;
        gr_ec_aff_point_t Pa, Aa;
        fmpz_t m, n, mn;
        slong j, nsmall, bits;
        int status;

        gr_ec_test_ring(R, state);

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        gr_ec_jac_point_init(P, E);
        gr_ec_jac_point_init(A, E);
        gr_ec_jac_point_init(B, E);
        gr_ec_jac_point_init(C, E);
        gr_ec_aff_point_init(Pa, E);
        gr_ec_aff_point_init(Aa, E);
        fmpz_init(m);
        fmpz_init(n);
        fmpz_init(mn);

        status = gr_ec_jac_point_randtest(P, state, E);

        if (status != GR_SUCCESS)
        {
            if (status & GR_UNABLE)
                count_unable++;
            else
                count_domain++;

            goto cleanup;
        }

        count_success++;

        bits = gr_ec_test_scalar_bits(R);
        fmpz_randtest(m, state, bits);
        fmpz_randtest(n, state, bits);
        fmpz_add(mn, m, n);

        /* 0 P = O and 1 P = P */
        GR_EC_TEST(gr_ec_jac_point_mul_ui(A, P, 0, E) == GR_SUCCESS, E);
        GR_EC_TEST(gr_ec_jac_point_is_inf(A, E) == T_TRUE, E);

        if (gr_ec_jac_point_mul_ui(A, P, 1, E) == GR_SUCCESS)
            GR_EC_TEST(gr_ec_jac_point_equal(A, P, E) != T_FALSE, E);

        /* a small multiple against repeated addition; the accumulator is
           never normalized, so its coordinates grow with every step and
           the chain has to stay short over an infinite ring */
        nsmall = 1 + n_randint(state,
                    (gr_ctx_is_finite(R) == T_TRUE) ? 12 : 4);

        if (gr_ec_jac_point_zero(B, E) == GR_SUCCESS)
        {
            status = GR_SUCCESS;

            for (j = 0; j < nsmall && status == GR_SUCCESS; j++)
                status = gr_ec_jac_point_add(B, B, P, E);

            if (status == GR_SUCCESS
                    && gr_ec_jac_point_mul_si(A, P, nsmall, E) == GR_SUCCESS)
            {
                GR_EC_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE, E);
                GR_EC_TEST(gr_ec_jac_point_is_on_curve(A, E) != T_FALSE, E);
            }
        }

        /* (m + n) P = m P + n P */
        if (gr_ec_jac_point_mul_fmpz(A, P, m, E) == GR_SUCCESS
                && gr_ec_jac_point_mul_fmpz(B, P, n, E) == GR_SUCCESS
                && gr_ec_jac_point_add(B, A, B, E) == GR_SUCCESS
                && gr_ec_jac_point_mul_fmpz(C, P, mn, E) == GR_SUCCESS)
            GR_EC_TEST(gr_ec_jac_point_equal(B, C, E) != T_FALSE, E);

        /* m (n P) = (m n) P */
        fmpz_mul(mn, m, n);

        if (gr_ec_jac_point_mul_fmpz(A, P, n, E) == GR_SUCCESS
                && gr_ec_jac_point_mul_fmpz(A, A, m, E) == GR_SUCCESS
                && gr_ec_jac_point_mul_fmpz(B, P, mn, E) == GR_SUCCESS)
            GR_EC_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE, E);

        /* (-n) P = -(n P) */
        fmpz_neg(mn, n);

        if (gr_ec_jac_point_mul_fmpz(A, P, mn, E) == GR_SUCCESS
                && gr_ec_jac_point_mul_fmpz(B, P, n, E) == GR_SUCCESS
                && gr_ec_jac_point_neg(B, B, E) == GR_SUCCESS)
            GR_EC_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE, E);

        /* the plain ladder agrees with the dispatched implementation */
        if (gr_ec_jac_point_mul_fmpz(A, P, n, E) == GR_SUCCESS
                && _gr_ec_jac_point_mul_fmpz_binary(B, P, n, E) == GR_SUCCESS)
            GR_EC_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE, E);

        /* aliasing */
        if (gr_ec_jac_point_mul_fmpz(A, P, n, E) == GR_SUCCESS)
        {
            GR_EC_TEST(gr_ec_jac_point_set(B, P, E) == GR_SUCCESS, E);
            GR_EC_TEST(gr_ec_jac_point_mul_fmpz(B, B, n, E) == GR_SUCCESS, E);
            GR_EC_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE, E);
        }

        /* the affine functions agree */
        if (gr_ec_jac_point_mul_fmpz(A, P, n, E) == GR_SUCCESS
                && gr_ec_aff_point_set_jac_point(Pa, P, E) == GR_SUCCESS)
        {
            if (gr_ec_aff_point_mul_fmpz(Aa, Pa, n, E) == GR_SUCCESS)
            {
                GR_EC_TEST(gr_ec_jac_point_set_aff_point(B, Aa, E) == GR_SUCCESS, E);
                GR_EC_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE, E);
            }

            /* the affine ladder does an inversion per bit; only run it
               where the coordinates stay bounded */
            if (gr_ctx_is_finite(R) == T_TRUE
                    && _gr_ec_aff_point_mul_fmpz_binary(Aa, Pa, n, E) == GR_SUCCESS)
            {
                GR_EC_TEST(gr_ec_jac_point_set_aff_point(B, Aa, E) == GR_SUCCESS, E);
                GR_EC_TEST(gr_ec_jac_point_equal(A, B, E) != T_FALSE, E);
            }
        }

cleanup:
        fmpz_clear(m);
        fmpz_clear(n);
        fmpz_clear(mn);
        gr_ec_aff_point_clear(Pa, E);
        gr_ec_aff_point_clear(Aa, E);
        gr_ec_jac_point_clear(P, E);
        gr_ec_jac_point_clear(A, E);
        gr_ec_jac_point_clear(B, E);
        gr_ec_jac_point_clear(C, E);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
