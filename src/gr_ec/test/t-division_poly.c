/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_poly.h"
#include "t-helpers.h"

/*
    In the short model the implementation's b-invariant formulas must come
    out as the classical a-invariant ones:

        psi_2^2 = 4 x^3 + 4 a4 x + 4 a6
        Psi_3   = 3 x^4 + 6 a4 x^2 + 12 a6 x - a4^2
        Psi_4   = 2 (x^6 + 5 a4 x^4 + 20 a6 x^3 - 5 a4^2 x^2
                        - 4 a4 a6 x - a4^3 - 8 a6^2)
*/
static void
check_short_model_base_cases(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_poly_t got, want;
        gr_ptr a4, a6, t, u;
        int status = GR_SUCCESS;

        gr_ec_test_ring(R, state);

        GR_TMP_INIT4(a4, a6, t, u, R);

        if (gr_randtest(a4, state, R) != GR_SUCCESS
                || gr_randtest(a6, state, R) != GR_SUCCESS
                || gr_ec_ctx_init_short_weierstrass(E, R, a4, a6) != GR_SUCCESS)
        {
            GR_TMP_CLEAR4(a4, a6, t, u, R);
            gr_ctx_clear(R);
            continue;
        }

        gr_poly_init(got, R);
        gr_poly_init(want, R);

        /* psi_2^2 */
        status |= gr_ec_ctx_psi2_sqr(got, E);
        status |= gr_poly_zero(want, R);
        status |= gr_poly_set_coeff_si(want, 3, 4, R);
        status |= gr_mul_ui(t, a4, 4, R);
        status |= gr_poly_set_coeff_scalar(want, 1, t, R);
        status |= gr_mul_ui(t, a6, 4, R);
        status |= gr_poly_set_coeff_scalar(want, 0, t, R);

        if (status == GR_SUCCESS)
            FLINT_TEST(gr_poly_equal(got, want, R) != T_FALSE);

        /* Psi_3 */
        status |= gr_ec_ctx_division_poly(got, 3, E);
        status |= gr_poly_zero(want, R);
        status |= gr_poly_set_coeff_si(want, 4, 3, R);
        status |= gr_mul_ui(t, a4, 6, R);
        status |= gr_poly_set_coeff_scalar(want, 2, t, R);
        status |= gr_mul_ui(t, a6, 12, R);
        status |= gr_poly_set_coeff_scalar(want, 1, t, R);
        status |= gr_sqr(t, a4, R);
        status |= gr_neg(t, t, R);
        status |= gr_poly_set_coeff_scalar(want, 0, t, R);

        if (status == GR_SUCCESS)
            FLINT_TEST(gr_poly_equal(got, want, R) != T_FALSE);

        /* Psi_4 */
        status |= gr_ec_ctx_division_poly(got, 4, E);
        status |= gr_poly_zero(want, R);
        status |= gr_poly_set_coeff_si(want, 6, 2, R);
        status |= gr_mul_ui(t, a4, 10, R);
        status |= gr_poly_set_coeff_scalar(want, 4, t, R);
        status |= gr_mul_ui(t, a6, 40, R);
        status |= gr_poly_set_coeff_scalar(want, 3, t, R);
        status |= gr_sqr(t, a4, R);
        status |= gr_mul_ui(t, t, 10, R);
        status |= gr_neg(t, t, R);
        status |= gr_poly_set_coeff_scalar(want, 2, t, R);
        status |= gr_mul(t, a4, a6, R);
        status |= gr_mul_ui(t, t, 8, R);
        status |= gr_neg(t, t, R);
        status |= gr_poly_set_coeff_scalar(want, 1, t, R);
        /* constant term -2 a4^3 - 16 a6^2 */
        status |= gr_sqr(t, a4, R);
        status |= gr_mul(t, t, a4, R);
        status |= gr_mul_two(t, t, R);
        status |= gr_sqr(u, a6, R);
        status |= gr_mul_ui(u, u, 16, R);
        status |= gr_add(t, t, u, R);
        status |= gr_neg(t, t, R);
        status |= gr_poly_set_coeff_scalar(want, 0, t, R);

        if (status == GR_SUCCESS)
            FLINT_TEST(gr_poly_equal(got, want, R) != T_FALSE);

        gr_poly_clear(got, R);
        gr_poly_clear(want, R);
        GR_TMP_CLEAR4(a4, a6, t, u, R);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    deg Psi_n = (n^2 - 1) / 2 with leading coefficient n for odd n, and
    (n^2 - 4) / 2 with leading coefficient n / 2 for even n. Only checked
    where the leading coefficient cannot vanish, that is over a base ring
    of characteristic 0 or of characteristic greater than n.
*/
static void
check_degree_and_leading(flint_rand_t state)
{
    slong which;

    for (which = 0; which < 3; which++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        ulong n, nmax = 24;
        gr_ptr t;

        switch (which)
        {
            case 0: gr_ctx_init_fmpz(R); break;
            case 1: gr_ctx_init_fmpq(R); break;
            default: GR_MUST_SUCCEED(gr_ctx_init_nmod(R, 1000003)); break;
        }

        FLINT_TEST(gr_ec_ctx_init_si(E, R, 1, 2, 3, 4, 5) == GR_SUCCESS);
        GR_TMP_INIT(t, R);

        for (n = 1; n <= nmax; n++)
        {
            gr_poly_t p;
            slong want_deg = (n % 2) ? (slong) (n * n - 1) / 2
                                    : (slong) (n * n - 4) / 2;
            gr_poly_init(p, R);

            FLINT_TEST(gr_ec_ctx_division_poly(p, n, E) == GR_SUCCESS);
            FLINT_TEST(gr_poly_length(p, R) == want_deg + 1);

            FLINT_TEST(gr_set_ui(t, (n % 2) ? n : n / 2, R) == GR_SUCCESS);
            FLINT_TEST(gr_equal(GR_ENTRY(p->coeffs, want_deg, R->sizeof_elem),
                        t, R) == T_TRUE);

            gr_poly_clear(p, R);
        }

        GR_TMP_CLEAR(t, R);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    The m-torsion sits inside the n-torsion whenever m divides n, so
    Psi_m divides Psi_n. Independent of how either is computed.
*/
static void
check_divisibility(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_poly_t pm, pn, q, r;
        ulong m, n;

        GR_MUST_SUCCEED(gr_ctx_init_nmod(R,
                    n_randprime(state, 10 + n_randint(state, 20), 1)));

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        m = 2 + n_randint(state, 6);
        n = m * (1 + n_randint(state, 4));

        gr_poly_init(pm, R);
        gr_poly_init(pn, R);
        gr_poly_init(q, R);
        gr_poly_init(r, R);

        if (gr_ec_ctx_division_poly(pm, m, E) == GR_SUCCESS
                && gr_ec_ctx_division_poly(pn, n, E) == GR_SUCCESS
                && gr_poly_length(pm, R) > 0
                && gr_poly_divrem(q, r, pn, pm, R) == GR_SUCCESS)
        {
            FLINT_TEST(gr_poly_is_zero(r, R) == T_TRUE);
        }

        gr_poly_clear(pm, R);
        gr_poly_clear(pn, R);
        gr_poly_clear(q, R);
        gr_poly_clear(r, R);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    The defining property, over a finite field so that both sides are
    decidable. For a point P other than O:

        Psi_n(x(P)) = 0                  ==>  n P = O
        n P = O and 2 P != O             ==>  Psi_n(x(P)) = 0

    The exception for even n is exactly the 2-torsion, where the psi_2 we
    divided out vanishes.
*/
static void
check_torsion(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, Q, D;
        gr_poly_t psi;
        gr_ptr x, y, v;
        ulong n;
        slong t;

        /* small fields, including characteristic 2 and 3 */
        switch (n_randint(state, 4))
        {
            case 0:
                gr_ctx_init_fq_nmod(R, 2, 2 + n_randint(state, 4), "a");
                break;
            case 1:
                gr_ctx_init_fq_nmod(R, 3, 1 + n_randint(state, 3), "a");
                break;
            default:
                GR_MUST_SUCCEED(gr_ctx_init_nmod(R,
                            n_randprime(state, 4 + n_randint(state, 6), 1)));
                break;
        }

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        n = 2 + n_randint(state, 11);

        gr_poly_init(psi, R);
        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);
        gr_ec_point_init(D, E);
        GR_TMP_INIT3(x, y, v, R);

        if (gr_ec_ctx_division_poly(psi, n, E) == GR_SUCCESS)
        {
            for (t = 0; t < 12; t++)
            {
                truth_t nP_inf, twoP_inf, psi_zero;

                if (gr_ec_point_randtest(P, state, E) != GR_SUCCESS
                        || gr_ec_point_is_inf(P, E) != T_FALSE
                        || gr_ec_point_get_affine(x, y, P, E) != GR_SUCCESS
                        || gr_poly_evaluate(v, psi, x, R) != GR_SUCCESS
                        || gr_ec_point_mul_ui(Q, P, n, E) != GR_SUCCESS
                        || gr_ec_point_dbl(D, P, E) != GR_SUCCESS)
                    continue;

                psi_zero = gr_is_zero(v, R);
                nP_inf = gr_ec_point_is_inf(Q, E);
                twoP_inf = gr_ec_point_is_inf(D, E);

                if (psi_zero == T_TRUE)
                    FLINT_TEST(nP_inf != T_FALSE);

                if (nP_inf == T_TRUE && !(n % 2 == 0 && twoP_inf != T_FALSE))
                    FLINT_TEST(psi_zero != T_FALSE);
            }
        }

        GR_TMP_CLEAR3(x, y, v, R);
        gr_ec_point_clear(D, E);
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(P, E);
        gr_poly_clear(psi, R);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

/*
    The table and the single-index ladder are different code paths, so
    they check each other's bookkeeping.
*/
static void
check_vec(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_poly_struct * tab;
        gr_poly_t p;
        slong len, k;

        gr_ec_test_ring(R, state);

        if (gr_ec_test_curve(E, R, state) != GR_SUCCESS)
        {
            gr_ctx_clear(R);
            continue;
        }

        /* over an infinite ring the coefficients grow fast */
        len = 1 + n_randint(state,
                    (gr_ctx_is_finite(R) == T_TRUE) ? 26 : 9);

        tab = flint_malloc(len * sizeof(gr_poly_struct));
        for (k = 0; k < len; k++)
            gr_poly_init(tab + k, R);

        gr_poly_init(p, R);

        if (gr_ec_ctx_division_poly_vec(tab, len, E) == GR_SUCCESS)
        {
            for (k = 0; k < len; k++)
            {
                FLINT_TEST(gr_ec_ctx_division_poly(p, (ulong) k, E) == GR_SUCCESS);
                FLINT_TEST(gr_poly_equal(p, tab + k, R) != T_FALSE);
            }
        }

        gr_poly_clear(p, R);
        for (k = 0; k < len; k++)
            gr_poly_clear(tab + k, R);
        flint_free(tab);

        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }
}

TEST_FUNCTION_START(gr_ec_division_poly, state)
{
    check_short_model_base_cases(state);
    check_degree_and_leading(state);
    check_divisibility(state);
    check_torsion(state);
    check_vec(state);

    TEST_FUNCTION_END(state);
}
