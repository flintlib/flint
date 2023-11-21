/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

static void
_arb_inv_si(arb_t res, slong a, slong prec)
{
    arb_set_si(res, a);
    arb_inv(res, res, prec);
}

/* Defined in t-platt_multieval.c and t-platt_multieval_threaded.c */
#ifndef _arb_div_si_si
#define _arb_div_si_si _arb_div_si_si
static void
_arb_div_si_si(arb_t res, slong a, slong b, slong prec)
{
    arb_set_si(res, a);
    arb_div_si(res, res, b, prec);
}
#endif

static void
_check_containment(const char *name, const arb_t x, const char *s)
{
    arb_t u;
    slong prec = 300;

    arb_init(u);
    arb_set_str(u, s, prec);

    if (!arb_contains(u, x))
    {
        flint_printf("FAIL: %s\n\n", name);
        flint_printf("observed = "); arb_printn(x, 30, 0); flint_printf("\n\n");
        flint_printf("expected = "); arb_printn(u, 30, 0); flint_printf("\n\n");
        flint_abort();
    }

    arb_clear(u);
}

TEST_FUNCTION_START(acb_dirichlet_platt_multieval, state)
{
    slong iter;

    /* Check a specific combination of parameter values that is relatively fast
     * to evaluate and that has relatively tight bounds. */
    {
        slong A = 8;
        slong B = 128;
        slong N = A*B;
        slong K = 30;
        slong sigma = 63;
        slong prec = 128;
        fmpz_t J, T;
        arb_t h;
        arb_ptr vec;

        arb_init(h);
        fmpz_init(J);
        fmpz_init(T);

        fmpz_set_si(J, 1000);
        fmpz_set_si(T, 10000);
        arb_set_d(h, 4.5);

        /* Spot-check lemma bound containment
         * in intervals calculated with PARI/GP. */
        {
            arb_t lem, xi, x, beta, t0;
            slong i = 201;
            slong k = 5;
            slong wp = 300;

            arb_init(lem);
            arb_init(xi);
            arb_init(x);
            arb_init(t0);
            arb_init(beta);

            _arb_inv_si(xi, B, wp);
            arb_mul_2exp_si(xi, xi, -1);
            _arb_div_si_si(x, i, B, wp);
            arb_set_fmpz(t0, T);
            acb_dirichlet_platt_beta(beta, t0, wp);

            acb_dirichlet_platt_lemma_32(lem, h, t0, x, wp);
            _check_containment("Lemma 3.2", lem,
                    "[5.3526496753240991744e-1072334 +/- 2.55e-1072354]");

            acb_dirichlet_platt_c_bound(lem, sigma, t0, h, k, wp);
            _check_containment("Lemma A.3", lem,
                    "[1.3516642396389823078e+134 +/- 2.65e+114]");

            acb_dirichlet_platt_lemma_A5(lem, B, h, k, wp);
            _check_containment("Lemma A.5", lem,
                    "[1.0075390047893384632e-30 +/- 5.57e-51]");

            acb_dirichlet_platt_lemma_A7(lem, sigma, t0, h, k, A, wp);
            _check_containment("Lemma A.7", lem,
                    "[3.0406705491484062400e-505 +/- 1.57e-525]");

            acb_dirichlet_platt_lemma_A9(lem, sigma, t0, h, A, wp);
            _check_containment("Lemma A.9", lem,
                    "[6.8953211848420326275e-536 +/- 3.52e-556]");

            acb_dirichlet_platt_lemma_A11(lem, t0, h, B, wp);
            _check_containment("Lemma A.11", lem,
                    "[3.0825745863006335768e-42 +/- 3.68e-62]");

            acb_dirichlet_platt_lemma_B1(lem, sigma, t0, h, J, wp);
            _check_containment("Lemma B.1", lem,
                    "[8.5737638613320328274e-42 +/- 7.50e-63]");

            acb_dirichlet_platt_lemma_B2(lem, K, h, xi, wp);
            _check_containment("Lemma B.2", lem,
                    "[2.0748437544358592615e-44 +/- 4.76e-64]");

            arb_clear(lem);
            arb_clear(xi);
            arb_clear(x);
            arb_clear(t0);
            arb_clear(beta);
        }

        /* Check a few random entries in the multieval vector. */
        vec = _arb_vec_init(N);
        acb_dirichlet_platt_multieval(vec, T, A, B, h, J, K, sigma, prec);
        for (iter = 0; iter < 20; iter++)
        {
            arb_t t, r;
            slong i = n_randint(state, N);
            slong n = i - N/2;
            arb_init(t);
            arb_init(r);
            _arb_div_si_si(t, n, A, prec);
            arb_add_fmpz(t, t, T, prec);
            acb_dirichlet_platt_scaled_lambda(r, t, prec);
            if (!arb_overlaps(vec + i, r))
            {
                flint_printf("FAIL: overlap for hardcoded example\n\n");
                flint_printf("i = %wd  n = %wd\n\n", i, n);
                flint_printf("vec[%wd] = ", i); arb_printn(vec + i, 30, 0); flint_printf("\n\n");
                flint_printf("r = "); arb_printn(r, 30, 0); flint_printf("\n\n");
                flint_abort();
            }
            arb_clear(t);
            arb_clear(r);
        }

        fmpz_clear(J);
        fmpz_clear(T);
        arb_clear(h);
        _arb_vec_clear(vec, N);
    }

    for (iter = 0; iter < 10 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec;
        ulong A, B, N, K;
        slong sigma, Tbits;
        fmpz_t J, T;
        arb_t h;
        arb_ptr v1, v2;

        /* better but slower limits are in parentheses below */
        prec = 2 + n_randint(state, 300);
        sigma = 1 + 2*(1 + n_randint(state, 100)); /* (200) */
        K = 1 + n_randint(state, 20); /* (50) */
        A = 1 + n_randint(state, 10);
        B = 1 + n_randint(state, 10); /* (500) */
        if (n_randint(state, 2))
            A *= 2;
        else
            B *= 2;
        N = A*B;

        fmpz_init(J);
        fmpz_init(T);
        fmpz_set_si(J, 1 + n_randint(state, 100)); /* (10000) */
        Tbits = 5 + n_randint(state, 15);
        fmpz_set_ui(T, n_randtest_bits(state, Tbits));

        arb_init(h);
        arb_set_si(h, 1 + n_randint(state, 20000));
        arb_div_si(h, h, 1000, prec);

        v1 = _arb_vec_init(N);
        v2 = _arb_vec_init(N);
        acb_dirichlet_platt_scaled_lambda_vec(v1, T, A, B, prec);
        acb_dirichlet_platt_multieval(v2, T, A, B, h, J, K, sigma, prec);

        if (!_arb_vec_overlaps(v1, v2, N))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd  prec = %wd\n\n", iter, prec);
            flint_printf("sigma = %wd\n\n", sigma);
            flint_printf("A = %wu  B = %wu  K = %wu\n\n", A, B, K);
            flint_printf("J = "); fmpz_print(J); flint_printf("\n\n");
            flint_printf("T = "); fmpz_print(T); flint_printf("\n\n");
            flint_printf("h = "); arb_printn(h, 30, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(h);
        fmpz_clear(J);
        fmpz_clear(T);
        _arb_vec_clear(v1, N);
        _arb_vec_clear(v2, N);
    }

    TEST_FUNCTION_END(state);
}
