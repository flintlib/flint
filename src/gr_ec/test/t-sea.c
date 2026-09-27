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
#include "fmpz_poly.h"
#include "gr_poly.h"
#include "gr_ec.h"
#include "gr_ec/impl.h"

static void
_init_phi(gr_poly_struct ** phi, ulong l, gr_ctx_t R)
{
    slong k;

    *phi = flint_malloc((l + 2) * sizeof(gr_poly_struct));

    for (k = 0; k < (slong) l + 2; k++)
        gr_poly_init(*phi + k, R);
}

static void
_clear_phi(gr_poly_struct * phi, ulong l, gr_ctx_t R)
{
    slong k;

    for (k = 0; k < (slong) l + 2; k++)
        gr_poly_clear(phi + k, R);

    flint_free(phi);
}

/* phi[a] == the integer polynomial given by its coefficients */
static int
_coeffs_are(const gr_poly_t f, const char ** c, slong n, gr_ctx_t R)
{
    fmpz_poly_t z;
    gr_poly_t g;
    fmpz_t x;
    slong i;
    int ok;

    fmpz_poly_init(z);
    fmpz_init(x);
    gr_poly_init(g, R);

    for (i = 0; i < n; i++)
    {
        fmpz_set_str(x, c[i], 10);
        fmpz_poly_set_coeff_fmpz(z, i, x);
    }

    GR_MUST_SUCCEED(gr_poly_set_fmpz_poly(g, z, R));
    ok = (gr_poly_equal(f, g, R) == T_TRUE);

    gr_poly_clear(g, R);
    fmpz_clear(x);
    fmpz_poly_clear(z);

    return ok;
}

/* Phi_2 and Phi^c_3, Phi^c_5 against their known integer coefficients */
static void
check_known(void)
{
    gr_ctx_t R;
    gr_poly_struct * phi;
    fmpz_t p;
    ulong s;

    fmpz_init(p);
    fmpz_set_str(p, "170141183460469231731687303715884105727", 10);
    gr_ctx_init_fmpz_mod(R, p);

    {
        const char * c3[] = {"1"};
        const char * c2[] = {"-162000", "1488", "-1"};
        const char * c1[] = {"8748000000", "40773375", "1488"};
        const char * c0[] = {"-157464000000000", "8748000000", "-162000", "1"};

        _init_phi(&phi, 2, R);
        FLINT_TEST(_gr_ec_modular_polynomial(phi, 2, R) == GR_SUCCESS);
        FLINT_TEST(_coeffs_are(phi + 3, c3, 1, R));
        FLINT_TEST(_coeffs_are(phi + 2, c2, 3, R));
        FLINT_TEST(_coeffs_are(phi + 1, c1, 3, R));
        FLINT_TEST(_coeffs_are(phi + 0, c0, 4, R));
        _clear_phi(phi, 2, R);
    }

    {
        const char * one[] = {"1"};
        const char * c3[] = {"36"}, * c2[] = {"270"};
        const char * c1[] = {"756", "-1"}, * c0[] = {"729"};

        _init_phi(&phi, 3, R);
        FLINT_TEST(_gr_ec_modular_polynomial_canonical(phi, &s, 3, R) == GR_SUCCESS);
        FLINT_TEST(s == 6);
        FLINT_TEST(_coeffs_are(phi + 4, one, 1, R));
        FLINT_TEST(_coeffs_are(phi + 3, c3, 1, R));
        FLINT_TEST(_coeffs_are(phi + 2, c2, 1, R));
        FLINT_TEST(_coeffs_are(phi + 1, c1, 2, R));
        FLINT_TEST(_coeffs_are(phi + 0, c0, 1, R));
        _clear_phi(phi, 3, R);
    }

    {
        const char * one[] = {"1"};
        const char * c5[] = {"30"}, * c4[] = {"315"}, * c3[] = {"1300"};
        const char * c2[] = {"1575"}, * c1[] = {"750", "-1"}, * c0[] = {"125"};

        _init_phi(&phi, 5, R);
        FLINT_TEST(_gr_ec_modular_polynomial_canonical(phi, &s, 5, R) == GR_SUCCESS);
        FLINT_TEST(s == 3);
        FLINT_TEST(_coeffs_are(phi + 6, one, 1, R));
        FLINT_TEST(_coeffs_are(phi + 5, c5, 1, R));
        FLINT_TEST(_coeffs_are(phi + 4, c4, 1, R));
        FLINT_TEST(_coeffs_are(phi + 3, c3, 1, R));
        FLINT_TEST(_coeffs_are(phi + 2, c2, 1, R));
        FLINT_TEST(_coeffs_are(phi + 1, c1, 2, R));
        FLINT_TEST(_coeffs_are(phi + 0, c0, 1, R));
        _clear_phi(phi, 5, R);
    }

    gr_ctx_clear(R);
    fmpz_clear(p);
}

/*
    Phi^c_l(f(q), j(q)) = 0 as power series, with f = l^s q^v eta(q^l)^2s /
    eta(q)^2s: multiplied by q^v, sum_a f^a sum_d c_ad q^(v - d) (q j)^d.
*/
static void
check_canonical_identity(flint_rand_t state)
{
    ulong ls[] = {7, 11, 13, 17, 19, 23, 29};
    slong li;

    for (li = 0; li < 7; li++)
    {
        ulong l = ls[li], s;
        slong v, N = 3 * l + 10, a, d;
        gr_ctx_t R;
        gr_poly_struct * phi;
        gr_poly_t f, J, fa, Jd, acc, term, num, den;
        fmpz_poly_t z;
        gr_ptr c;

        while (gr_ctx_init_nmod(R, n_randprime(state, 30 + n_randint(state, 30), 1)) != GR_SUCCESS)
            ;

        _init_phi(&phi, l, R);
        FLINT_TEST(_gr_ec_modular_polynomial_canonical(phi, &s, l, R) == GR_SUCCESS);
        v = s * (l - 1) / 12;

        fmpz_poly_init(z);
        gr_poly_init(f, R); gr_poly_init(J, R); gr_poly_init(fa, R);
        gr_poly_init(Jd, R); gr_poly_init(acc, R); gr_poly_init(term, R);
        gr_poly_init(num, R); gr_poly_init(den, R);
        GR_TMP_INIT(c, R);

        fmpz_poly_eta_qexp(z, 2 * s, N);
        GR_MUST_SUCCEED(gr_poly_set_fmpz_poly(den, z, R));
        GR_MUST_SUCCEED(gr_poly_inv_series(den, den, N, R));
        fmpz_poly_eta_qexp(z, 2 * s, N / l + 1);
        fmpz_poly_inflate(z, z, l);
        GR_MUST_SUCCEED(gr_poly_set_fmpz_poly(num, z, R));
        GR_MUST_SUCCEED(gr_poly_mullow(f, num, den, N, R));
        GR_MUST_SUCCEED(gr_set_ui(c, l, R));
        GR_MUST_SUCCEED(gr_pow_ui(c, c, s, R));
        GR_MUST_SUCCEED(gr_poly_mul_scalar(f, f, c, R));
        GR_MUST_SUCCEED(gr_poly_shift_left(f, f, v, R));
        GR_MUST_SUCCEED(_gr_ec_j_qexp(J, N, R));

        GR_MUST_SUCCEED(gr_poly_zero(acc, R));
        GR_MUST_SUCCEED(gr_poly_one(fa, R));

        for (a = 0; a < (slong) l + 2; a++)
        {
            GR_MUST_SUCCEED(gr_poly_one(Jd, R));

            for (d = 0; d < phi[a].length; d++)
            {
                GR_MUST_SUCCEED(gr_poly_mul_scalar(term, Jd, GR_ENTRY(phi[a].coeffs, d, R->sizeof_elem), R));
                GR_MUST_SUCCEED(gr_poly_shift_left(term, term, v - d, R));
                GR_MUST_SUCCEED(gr_poly_mullow(term, term, fa, N, R));
                GR_MUST_SUCCEED(gr_poly_add(acc, acc, term, R));
                GR_MUST_SUCCEED(gr_poly_mullow(Jd, Jd, J, N, R));
            }

            GR_MUST_SUCCEED(gr_poly_mullow(fa, fa, f, N, R));
        }

        GR_MUST_SUCCEED(gr_poly_truncate(acc, acc, N - v - 1, R));
        FLINT_TEST(gr_poly_is_zero(acc, R) == T_TRUE);

        GR_TMP_CLEAR(c, R);
        gr_poly_clear(f, R); gr_poly_clear(J, R); gr_poly_clear(fa, R);
        gr_poly_clear(Jd, R); gr_poly_clear(acc, R); gr_poly_clear(term, R);
        gr_poly_clear(num, R); gr_poly_clear(den, R);
        fmpz_poly_clear(z);
        _clear_phi(phi, l, R);
        gr_ctx_clear(R);
    }
}

static int
_random_curve(gr_ec_ctx_t E, gr_ctx_t R, flint_rand_t state)
{
    gr_ptr a, b;
    slong tries;
    int ok = 0;

    GR_TMP_INIT2(a, b, R);

    for (tries = 0; tries < 50 && !ok; tries++)
        if (gr_randtest_not_zero(a, state, R) == GR_SUCCESS
                && gr_randtest_not_zero(b, state, R) == GR_SUCCESS
                && gr_ec_ctx_init_short_weierstrass(E, R, a, b) == GR_SUCCESS)
            ok = 1;

    GR_TMP_CLEAR2(a, b, R);

    return ok;
}

/*
    At every prime l: the Elkies trace equals t mod l, and an Atkin set
    contains t mod l.
*/
static void
check_elkies_atkin(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t p, N, t;
        ulong l;

        fmpz_init(p); fmpz_init(N); fmpz_init(t);
        fmpz_randprime(p, state, 30 + n_randint(state, 20), 0);
        gr_ctx_init_fmpz_mod(R, p);
        GR_MUST_SUCCEED(gr_ctx_set_is_field(R, T_TRUE));

        if (_random_curve(E, R, state))
        {
            FLINT_TEST(gr_ec_ctx_cardinality_bsgs(N, E) == GR_SUCCESS);
            fmpz_add_ui(t, p, 1);
            fmpz_sub(t, t, N);

            for (l = 3; l < 40; l = n_nextprime(l, 1))
            {
                ulong tl, r, truth = fmpz_fdiv_ui(t, l), T[64];
                int el;
                slong n, k, in = 0;

                if (_gr_ec_elkies_trace(&tl, &el, &r, l, p, E) != GR_SUCCESS)
                    continue;

                if (el)
                    FLINT_TEST(tl == truth);
                else if (r != 0)
                {
                    n = _gr_ec_atkin_candidates(T, l, r, fmpz_fdiv_ui(p, l));

                    for (k = 0; k < n; k++)
                        in |= (T[k] == truth);

                    FLINT_TEST(in);
                }
            }

            gr_ec_ctx_clear(E);
        }

        gr_ctx_clear(R);
        fmpz_clear(p); fmpz_clear(N); fmpz_clear(t);
    }
}

/*
    Match and sort on its own: t mod m3 and made up candidate sets that
    contain the truth, with decoys, must give back the order.
*/
static void
check_match_sort(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t p, N, N2, t, m3, t3;
        gr_ec_atkin_struct A[4];
        ulong ls[4] = {13, 17, 19, 23};
        slong i, k;

        fmpz_init(p); fmpz_init(N); fmpz_init(N2); fmpz_init(t);
        fmpz_init(m3); fmpz_init(t3);
        fmpz_randprime(p, state, 24 + n_randint(state, 24), 0);
        gr_ctx_init_fmpz_mod(R, p);
        GR_MUST_SUCCEED(gr_ctx_set_is_field(R, T_TRUE));

        if (_random_curve(E, R, state))
        {
            FLINT_TEST(gr_ec_ctx_cardinality_bsgs(N, E) == GR_SUCCESS);
            fmpz_add_ui(t, p, 1);
            fmpz_sub(t, t, N);

            fmpz_set_ui(m3, 1 + n_randint(state, 60));
            fmpz_mod(t3, t, m3);

            for (i = 0; i < 4; i++)
            {
                ulong l = ls[i], truth = fmpz_fdiv_ui(t, l);

                A[i].l = l;
                A[i].T = flint_malloc(l * sizeof(ulong));
                A[i].n = 1 + n_randint(state, 3);
                A[i].T[0] = truth;

                for (k = 1; k < A[i].n; k++)
                    A[i].T[k] = (truth + k) % l;
            }

            FLINT_TEST(_gr_ec_cardinality_match_sort(N2, t3, m3, A, 4, E) == GR_SUCCESS);
            FLINT_TEST(fmpz_equal(N, N2));

            for (i = 0; i < 4; i++)
                flint_free(A[i].T);

            gr_ec_ctx_clear(E);
        }

        gr_ctx_clear(R);
        fmpz_clear(p); fmpz_clear(N); fmpz_clear(N2); fmpz_clear(t);
        fmpz_clear(m3); fmpz_clear(t3);
    }
}

/* SEA against BSGS and the naive count, over prime and extension fields */
static void
check_sea(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t N1, N2, N3, p;
        int kind = n_randint(state, 3);

        fmpz_init(N1); fmpz_init(N2); fmpz_init(N3); fmpz_init(p);

        if (kind == 0)
        {
            while (gr_ctx_init_nmod(R, n_randprime(state, 8 + n_randint(state, 12), 1)) != GR_SUCCESS)
                ;
        }
        else if (kind == 1)
        {
            fmpz_randprime(p, state, 30 + n_randint(state, 20), 0);
            gr_ctx_init_fmpz_mod(R, p);
            GR_MUST_SUCCEED(gr_ctx_set_is_field(R, T_TRUE));
        }
        else
        {
            gr_ctx_init_fq_nmod(R, n_randprime(state, 5 + n_randint(state, 5), 1),
                    2 + n_randint(state, 2), "a");
        }

        if (_random_curve(E, R, state))
        {
            if (gr_ec_ctx_cardinality_sea(N1, E) == GR_SUCCESS)
            {
                if (gr_ec_ctx_cardinality_bsgs(N2, E) == GR_SUCCESS)
                    FLINT_TEST(fmpz_equal(N1, N2));

                if (gr_ctx_cardinality_fmpz(p, R) == GR_SUCCESS
                        && fmpz_cmp_ui(p, 70000) < 0
                        && gr_ec_ctx_cardinality_naive(N3, E) == GR_SUCCESS)
                    FLINT_TEST(fmpz_equal(N1, N3));
            }

            gr_ec_ctx_clear(E);
        }

        gr_ctx_clear(R);
        fmpz_clear(N1); fmpz_clear(N2); fmpz_clear(N3); fmpz_clear(p);
    }
}

TEST_FUNCTION_START(gr_ec_sea, state)
{
    check_known();
    check_canonical_identity(state);
    check_elkies_atkin(state);
    check_match_sort(state);
    check_sea(state);

    TEST_FUNCTION_END(state);
}
