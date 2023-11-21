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

TEST_FUNCTION_START(acb_dirichlet_platt_ws_interpolation, state)
{
    slong iter;

    for (iter = 0; iter < 5 * 0.1 * flint_test_multiplier(); iter++)
    {
        ulong A, B, iter2;
        int Tbits;
        fmpz_t T;
        arb_ptr vec;
        slong vprec;

        A = 1 + n_randint(state, 10);
        B = 1 + n_randint(state, 10);
        if (n_randint(state, 2))
            A *= 2;
        else
            B *= 2;
        fmpz_init(T);
        Tbits = 5 + n_randint(state, 20);
        fmpz_set_ui(T, n_randtest_bits(state, Tbits));
        vprec = 2 + n_randint(state, 300);

        vec = _arb_vec_init(A*B);
        acb_dirichlet_platt_scaled_lambda_vec(vec, T, A, B, vprec);

        for (iter2 = 0; iter2 < 50; iter2++)
        {
            ulong Ns_max, sigma;
            arb_t expected, observed;
            arb_t t0, H;
            slong prec;

            arb_init(t0);
            arb_init(H);
            arb_init(expected);
            arb_init(observed);

            prec = 2 + n_randint(state, 300);
            Ns_max = 1 + n_randint(state, 100);
            sigma = 1 + 2*(1 + n_randint(state, 100));
            arb_set_si(t0, A*B*500 - n_randint(state, A*B*1000) - 1);
            arb_div_ui(t0, t0, A*1000, prec);
            arb_add_fmpz(t0, t0, T, prec);
            arb_set_ui(H, 1 + n_randint(state, 10000));
            arb_div_ui(H, H, 1000, prec);
            arb_abs(H, H);

            acb_dirichlet_platt_scaled_lambda(expected, t0, prec);
            acb_dirichlet_platt_ws_interpolation(observed, NULL, t0, vec,
                    T, A, B, Ns_max, H, sigma, prec);

            if (!arb_overlaps(expected, observed))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd\n\n", iter);
                flint_printf("iter2 = %wd\n\n", iter2);
                flint_printf("A = %wu\n\n", A);
                flint_printf("B = %wu\n\n", B);
                flint_printf("T = "); fmpz_print(T); flint_printf("\n\n");
                flint_printf("vprec = %wd\n\n", vprec);
                flint_printf("Ns_max = %wu\n\n", Ns_max);
                flint_printf("sigma = %wu\n\n", sigma);
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("t0 = "); arb_printd(t0, 15); flint_printf("\n\n");
                flint_printf("H = "); arb_printd(H, 15); flint_printf("\n\n");
            }

            arb_clear(t0);
            arb_clear(H);
            arb_clear(expected);
            arb_clear(observed);
        }

        fmpz_clear(T);
        _arb_vec_clear(vec, A*B);
    }

    TEST_FUNCTION_END(state);
}
