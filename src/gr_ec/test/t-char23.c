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
    Curves over fields of characteristic 2 and 3. In characteristic 2 the
    short Weierstrass equation is always singular, so the general form is
    the only one that represents a curve at all, and the long Weierstrass
    formulas are not merely a generalization of the short ones but the
    only thing that works.
*/
TEST_GR_FUNCTION_START(gr_ec_char23, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, Q, S, A, B;
        gr_ec_jac_point_t Pj, Qj, Sj;
        gr_ec_aff_point_t Pa, Qa, Sa;
        fmpz_t m, n;
        int status;

        gr_ctx_init_fq_nmod(R, 2 + n_randint(state, 2), 1 + n_randint(state, 5), "a");

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        FLINT_TEST(gr_ec_ctx_is_smooth(E) == T_TRUE);

        /* In characteristic 2 a short Weierstrass equation is always
           singular, so the context can only have picked the long model.
           In characteristic 3 it is nonsingular whenever a4 is nonzero,
           and both models occur. */
        if (gr_ctx_is_finite_characteristic(R) == T_TRUE
                && gr_ec_ctx_model(E) == GR_EC_SHORT_WEIERSTRASS)
            FLINT_TEST(gr_is_zero(GR_EC_A4(E), R) == T_FALSE);

        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);
        gr_ec_point_init(S, E);
        gr_ec_point_init(A, E);
        gr_ec_point_init(B, E);
        gr_ec_jac_point_init(Pj, E);
        gr_ec_jac_point_init(Qj, E);
        gr_ec_jac_point_init(Sj, E);
        gr_ec_aff_point_init(Pa, E);
        gr_ec_aff_point_init(Qa, E);
        gr_ec_aff_point_init(Sa, E);
        fmpz_init(m);
        fmpz_init(n);

        status = gr_ec_point_randtest(P, state, E);
        status |= gr_ec_point_randtest(Q, state, E);

        if (status != GR_SUCCESS)
        {
            count_unable++;
            goto cleanup;
        }

        count_success++;

        FLINT_TEST(gr_ec_point_is_on_curve(P, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_is_on_curve(Q, E) == T_TRUE);

        /* P - P = O */
        FLINT_TEST(gr_ec_point_neg(A, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_is_on_curve(A, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_add(B, P, A, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_is_inf(B, E) == T_TRUE);

        /* P + P = 2P */
        FLINT_TEST(gr_ec_point_add(A, P, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_dbl(B, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_is_on_curve(B, E) == T_TRUE);

        /* commutativity and associativity */
        FLINT_TEST(gr_ec_point_add(A, P, Q, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_add(B, Q, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_is_on_curve(A, E) == T_TRUE);

        FLINT_TEST(gr_ec_point_add(S, P, Q, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_add(A, S, Q, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_dbl(S, Q, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_add(B, P, S, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);

        /* the three representations agree */
        FLINT_TEST(gr_ec_point_add(S, P, Q, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_jac_point_set_point(Pj, P, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_jac_point_set_point(Qj, Q, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_jac_point_is_on_curve(Pj, E) == T_TRUE);
        FLINT_TEST(gr_ec_jac_point_add(Sj, Pj, Qj, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_set_jac_point(A, Sj, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(A, S, E) == T_TRUE);

        if (gr_ec_aff_point_set_point(Pa, P, E) == GR_SUCCESS
                && gr_ec_aff_point_set_point(Qa, Q, E) == GR_SUCCESS)
        {
            FLINT_TEST(gr_ec_aff_point_is_on_curve(Pa, E) == T_TRUE);

            FLINT_TEST(gr_ec_aff_point_add(Sa, Pa, Qa, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_set_aff_point(A, Sa, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_equal(A, S, E) == T_TRUE);

            FLINT_TEST(gr_ec_jac_point_add_aff_point(Sj, Pj, Qa, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_set_jac_point(A, Sj, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_point_equal(A, S, E) == T_TRUE);
        }

        /* scalar multiplication, which runs the NAF ladder */
        fmpz_randtest(m, state, 40);
        fmpz_randtest(n, state, 40);

        FLINT_TEST(gr_ec_point_mul_fmpz(A, P, m, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_mul_fmpz(B, P, n, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_add(S, A, B, E) == GR_SUCCESS);
        fmpz_add(m, m, n);
        FLINT_TEST(gr_ec_point_mul_fmpz(A, P, m, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(A, S, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_is_on_curve(A, E) == T_TRUE);

        /* and agrees with the plain ladder */
        FLINT_TEST(_gr_ec_point_mul_fmpz_binary(B, P, m, E) == GR_SUCCESS);
        FLINT_TEST(gr_ec_point_equal(A, B, E) == T_TRUE);

cleanup:
        fmpz_clear(m);
        fmpz_clear(n);
        gr_ec_aff_point_clear(Pa, E);
        gr_ec_aff_point_clear(Qa, E);
        gr_ec_aff_point_clear(Sa, E);
        gr_ec_jac_point_clear(Pj, E);
        gr_ec_jac_point_clear(Qj, E);
        gr_ec_jac_point_clear(Sj, E);
        gr_ec_point_clear(P, E);
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(S, E);
        gr_ec_point_clear(A, E);
        gr_ec_point_clear(B, E);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
