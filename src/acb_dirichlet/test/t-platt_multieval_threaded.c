/*
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

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

TEST_FUNCTION_START(acb_dirichlet_platt_multieval_threaded, state)
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
        fmpz_t J;
        fmpz_t T;
        arb_t h;
        arb_ptr vec;

        arb_init(h);
        fmpz_init(J);
        fmpz_init(T);

        fmpz_set_si(J, 1000);
        fmpz_set_si(T, 10000);
        arb_set_d(h, 4.5);

        flint_set_num_threads(5);

        /* Check a few random entries in the multieval vector. */
        vec = _arb_vec_init(N);
        acb_dirichlet_platt_multieval_threaded(vec, T, A, B, h, J, K, sigma, prec);
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

        flint_set_num_threads(1 + n_randint(state, 5));

        v1 = _arb_vec_init(N);
        v2 = _arb_vec_init(N);
        acb_dirichlet_platt_scaled_lambda_vec(v1, T, A, B, prec);
        acb_dirichlet_platt_multieval_threaded(v2, T, A, B, h, J, K, sigma, prec);

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
