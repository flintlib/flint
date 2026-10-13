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

#define CACHED(E) (&GR_EC_CTX(E)->cardinality)

TEST_FUNCTION_START(gr_ec_ctx_set_cardinality, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_jac_point_t P, A, B;
        fmpz_t N, M, k;
        ulong p;

        /* small fields are settled by counting, larger ones by probing */
        p = n_randprime(state, 5 + n_randint(state, 20), 1);

        if (p <= 3 || gr_ctx_init_nmod(R, p) != GR_SUCCESS)
            continue;

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(N);
        fmpz_init(M);
        fmpz_init(k);
        gr_ec_jac_point_init(P, E);
        gr_ec_jac_point_init(A, E);
        gr_ec_jac_point_init(B, E);

        FLINT_TEST(gr_ec_ctx_cardinality_bsgs(N, E) == GR_SUCCESS
            || gr_ec_ctx_cardinality_naive(N, E) == GR_SUCCESS);

        /* the generic cardinality agrees, and does not cache */
        FLINT_TEST(gr_ctx_cardinality_fmpz(M, E) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal(M, N));
        FLINT_TEST(fmpz_is_zero(CACHED(E)));

        /* a value outside the Hasse interval is refused */
        fmpz_mul_ui(M, N, 100);
        fmpz_add_ui(M, M, 1);
        FLINT_TEST(gr_ec_ctx_set_cardinality(E, M) == GR_DOMAIN);
        FLINT_TEST(fmpz_is_zero(CACHED(E)));

        fmpz_zero(M);
        FLINT_TEST(gr_ec_ctx_set_cardinality(E, M) == GR_DOMAIN);

        /* the truth is accepted and cached */
        FLINT_TEST(gr_ec_ctx_set_cardinality(E, N) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal(CACHED(E), N));
        FLINT_TEST(gr_ec_ctx_cardinality(M, E) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal(M, N));

        /* once known, anything else is refused */
        fmpz_add_ui(M, N, 1);
        FLINT_TEST(gr_ec_ctx_set_cardinality(E, M) == GR_DOMAIN);
        FLINT_TEST(fmpz_equal(CACHED(E), N));

        /* reducing scalars modulo the cached value changes nothing */
        if (gr_ec_jac_point_randtest(P, state, E) == GR_SUCCESS)
        {
            fmpz_randbits(k, state, 1 + n_randint(state, 200));
            FLINT_TEST(gr_ec_jac_point_mul_fmpz(A, P, k, E) == GR_SUCCESS);

            /* the same multiplication without the cache */
            fmpz_zero(CACHED(E));
            FLINT_TEST(gr_ec_jac_point_mul_fmpz(B, P, k, E) == GR_SUCCESS);
            FLINT_TEST(gr_ec_jac_point_equal(A, B, E) == T_TRUE);
        }

        /* gr_ec_ctx_cardinality counts and caches */
        fmpz_zero(CACHED(E));
        FLINT_TEST(gr_ec_ctx_cardinality(M, E) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal(M, N));
        FLINT_TEST(fmpz_equal(CACHED(E), N));

        /* a wrong value within the Hasse interval is caught by a point
           it does not kill, whenever there is one */
        fmpz_zero(CACHED(E));
        fmpz_add_ui(M, N, 1);

        if (gr_ec_jac_point_randtest(P, state, E) == GR_SUCCESS
                && gr_ec_jac_point_mul_fmpz(A, P, M, E) == GR_SUCCESS
                && gr_ec_jac_point_is_inf(A, E) == T_FALSE)
        {
            FLINT_TEST(gr_ec_ctx_set_cardinality(E, M) == GR_DOMAIN);
            FLINT_TEST(fmpz_is_zero(CACHED(E)));
        }

        gr_ec_jac_point_clear(P, E);
        gr_ec_jac_point_clear(A, E);
        gr_ec_jac_point_clear(B, E);
        fmpz_clear(N);
        fmpz_clear(M);
        fmpz_clear(k);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    /* E(R) is finite over a finite R, and not provably infinite otherwise:
       y^2 = x^3 - x has four points over Q */
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;

        gr_ctx_init_fmpq(R);
        FLINT_TEST(gr_ec_ctx_init_short_weierstrass_si(E, R, -1, 0) == GR_SUCCESS);
        FLINT_TEST(gr_ctx_is_finite(E) == T_UNKNOWN);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);

        FLINT_TEST(gr_ctx_init_nmod(R, 101) == GR_SUCCESS);
        FLINT_TEST(gr_ec_ctx_init_short_weierstrass_si(E, R, 1, 1) == GR_SUCCESS);
        FLINT_TEST(gr_ctx_is_finite(E) == T_TRUE);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}
