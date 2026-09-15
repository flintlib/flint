/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "t-helpers.h"

/* Sanity check on the generic test function itself, over groups whose
   additive structure is not in doubt. */
static void
check_additive_group_on_rings(flint_rand_t state)
{
    gr_ctx_t R;

    gr_ctx_init_fmpz(R);
    gr_test_additive_group(R, 100, GR_TEST_ALWAYS_ABLE);
    gr_ctx_clear(R);

    gr_ctx_init_fmpq(R);
    gr_test_additive_group(R, 100, GR_TEST_ALWAYS_ABLE);
    gr_ctx_clear(R);

    while (gr_ctx_init_nmod(R, n_randtest_not_zero(state)) != GR_SUCCESS)
        ;
    gr_test_additive_group(R, 100, GR_TEST_ALWAYS_ABLE);
    gr_ctx_clear(R);
}

/* The operations that do not exist on a curve must say so with GR_DOMAIN,
   not with GR_UNABLE: the caller has to be able to tell "there is no such
   thing here" from "I could not compute it". */
static void
check_ring_operations_are_out_of_domain(gr_ec_ctx_t E)
{
    gr_ptr P, Q;
    fmpz_t n;
    fmpq_t q;

    GR_TMP_INIT2(P, Q, E);
    fmpz_init(n);
    fmpq_init(q);

    FLINT_TEST(gr_one(P, E) == GR_DOMAIN);
    FLINT_TEST(gr_neg_one(P, E) == GR_DOMAIN);
    FLINT_TEST(gr_set_ui(P, 1, E) == GR_DOMAIN);
    FLINT_TEST(gr_set_si(P, 1, E) == GR_DOMAIN);
    FLINT_TEST(gr_set_fmpz(P, n, E) == GR_DOMAIN);
    FLINT_TEST(gr_set_fmpq(P, q, E) == GR_DOMAIN);
    FLINT_TEST(gr_mul(P, P, Q, E) == GR_DOMAIN);
    FLINT_TEST(gr_sqr(P, Q, E) == GR_DOMAIN);
    FLINT_TEST(gr_mul_fmpq(P, Q, q, E) == GR_DOMAIN);
    FLINT_TEST(gr_div(P, P, Q, E) == GR_DOMAIN);
    FLINT_TEST(gr_inv(P, Q, E) == GR_DOMAIN);
    FLINT_TEST(gr_pow_ui(P, Q, 2, E) == GR_DOMAIN);
    FLINT_TEST(gr_pow_si(P, Q, 2, E) == GR_DOMAIN);
    FLINT_TEST(gr_pow_fmpz(P, Q, n, E) == GR_DOMAIN);

    /* is_one has no answer rather than a wrong one */
    FLINT_TEST(gr_is_one(P, E) != T_TRUE);

    /* and the domain does not claim to be something it is not */
    FLINT_TEST(gr_ctx_is_ring(E) == T_FALSE);
    FLINT_TEST(gr_ctx_is_field(E) == T_FALSE);
    FLINT_TEST(gr_ctx_is_multiplicative_group(E) == T_FALSE);
    FLINT_TEST(gr_ctx_base(E) == (gr_ptr) GR_EC_ELEM_CTX(E));

    fmpq_clear(q);
    fmpz_clear(n);
    GR_TMP_CLEAR2(P, Q, E);
}

/* The generic interface and the module functions must be the same thing.
   Only the projective representation is checked directly; the other two
   are covered by t-convert plus the group tests. */
static void
check_generic_matches_module(gr_ec_ctx_t E, flint_rand_t state)
{
    gr_ec_point_t P, Q, res1, res2;
    fmpz_t n;
    int status;

    gr_ec_point_init(P, E);
    gr_ec_point_init(Q, E);
    gr_ec_point_init(res1, E);
    gr_ec_point_init(res2, E);
    fmpz_init(n);

    status = gr_ec_point_randtest(P, state, E);
    status |= gr_ec_point_randtest(Q, state, E);

    fmpz_randtest(n, state, gr_ec_test_scalar_bits(GR_EC_ELEM_CTX(E)));

    if (status == GR_SUCCESS)
    {
        status = gr_add(res1, P, Q, E);
        status |= gr_ec_point_add(res2, P, Q, E);
        if (status == GR_SUCCESS)
            FLINT_TEST(gr_ec_point_equal(res1, res2, E) != T_FALSE);

        status = gr_sub(res1, P, Q, E);
        status |= gr_ec_point_sub(res2, P, Q, E);
        if (status == GR_SUCCESS)
            FLINT_TEST(gr_ec_point_equal(res1, res2, E) != T_FALSE);

        status = gr_neg(res1, P, E);
        status |= gr_ec_point_neg(res2, P, E);
        if (status == GR_SUCCESS)
            FLINT_TEST(gr_ec_point_equal(res1, res2, E) != T_FALSE);

        status = gr_mul_fmpz(res1, P, n, E);
        status |= gr_ec_point_mul_fmpz(res2, P, n, E);
        if (status == GR_SUCCESS)
            FLINT_TEST(gr_ec_point_equal(res1, res2, E) != T_FALSE);

        FLINT_TEST(gr_zero(res1, E) == GR_SUCCESS);
        FLINT_TEST(gr_is_zero(res1, E) == T_TRUE);
        FLINT_TEST(gr_ec_point_is_inf(res1, E) == T_TRUE);
    }

    fmpz_clear(n);
    gr_ec_point_clear(res2, E);
    gr_ec_point_clear(res1, E);
    gr_ec_point_clear(Q, E);
    gr_ec_point_clear(P, E);
}

TEST_FUNCTION_START(gr_ec_generic, state)
{
    slong iter;

    check_additive_group_on_rings(state);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_repr_t repr;
        int flags = 0;

        gr_ec_test_ring(R, state);

        repr = n_randint(state, GR_EC_NUM_REPRS);

        /* affine arithmetic needs to invert */
        if (repr == GR_EC_REPR_AFFINE && gr_ctx_is_field(R) != T_TRUE)
            repr = GR_EC_REPR_PROJECTIVE;

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        FLINT_TEST(gr_ec_ctx_repr(E) == GR_EC_REPR_PROJECTIVE);
        FLINT_TEST(gr_ec_ctx_set_repr(E, repr) == GR_SUCCESS);
        FLINT_TEST(gr_ec_ctx_repr(E) == repr);

        check_ring_operations_are_out_of_domain(E);

        /* Over a finite field every operation is decidable, so nothing is
           allowed to come back GR_UNABLE there. */
        if (gr_ctx_is_finite(R) == T_TRUE && gr_ctx_is_field(R) == T_TRUE)
            flags = GR_TEST_ALWAYS_ABLE;

        gr_test_additive_group(E, 10, flags);

        if (repr == GR_EC_REPR_PROJECTIVE)
            check_generic_matches_module(E, state);

        /* gr_ctx_clear dispatches to gr_ec_ctx_clear */
        gr_ctx_clear(E);
        gr_ctx_clear(R);
    }

    /* the constructor that takes the representation up front */
    {
        gr_ctx_t R, E;

        gr_ctx_init_fmpq(R);

        FLINT_TEST(gr_ctx_init_gr_ec(E, R, NULL, NULL, NULL, NULL, NULL,
                    GR_EC_NUM_REPRS) == GR_DOMAIN);

        gr_ctx_clear(R);

        /* Z/4 is not a field, so affine is refused */
        while (gr_ctx_init_nmod(R, 4) != GR_SUCCESS)
            ;
        FLINT_TEST(gr_ctx_is_field(R) == T_FALSE);

        {
            gr_ptr a;
            GR_TMP_INIT(a, R);
            FLINT_TEST(gr_set_si(a, 1, R) == GR_SUCCESS);
            FLINT_TEST(gr_ctx_init_gr_ec(E, R, a, a, a, a, a,
                        GR_EC_REPR_AFFINE) == GR_DOMAIN);
            GR_TMP_CLEAR(a, R);
        }

        gr_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}
