/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

/* Conflicting types in t-rising_ui.c and t-rising_ui_jet.c */
#define rising_algorithm rising_algorithm1
void
rising_algorithm(acb_ptr res, const acb_t x, ulong n, ulong m, slong len, slong prec, int alg)
{
    if (alg == 0)
        acb_hypgeom_rising_ui_jet_powsum(res, x, n, len, prec);
    else if (alg == 1)
        acb_hypgeom_rising_ui_jet_rs(res, x, n, m, len, prec);
    else if (alg == 2)
        acb_hypgeom_rising_ui_jet_bs(res, x, n, len, prec);
    else
        acb_hypgeom_rising_ui_jet(res, x, n, len, prec);
}

TEST_FUNCTION_START(acb_hypgeom_rising_ui_jet, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, xk;
        acb_ptr y, ya, yb, yayb;
        ulong k, n, m1, m2, m3, len;
        slong prec;
        int alg1, alg2, alg3;

        prec = 2 + n_randint(state, 200);
        len = 1 + n_randint(state, 6);
        k = n_randint(state, 10);
        n = n_randint(state, 50);
        m1 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n + k, 1));
        m2 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(k, 1));
        m3 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n, 1));
        alg1 = n_randint(state, 4);
        alg2 = n_randint(state, 4);
        alg3 = n_randint(state, 4);

        acb_init(x);
        acb_init(xk);
        y = _acb_vec_init(len);
        ya = _acb_vec_init(len);
        yb = _acb_vec_init(len);
        yayb = _acb_vec_init(len);

        acb_randtest(x, state, prec, 10);
        acb_add_ui(xk, x, k, prec);

        rising_algorithm(y, x, n + k, m1, len, prec, alg1);
        rising_algorithm(ya, x, k, m2, len, prec, alg2);
        rising_algorithm(yb, xk, n, m3, len, prec, alg3);
        _acb_poly_mullow(yayb, ya, len, yb, len, len, prec);

        if (!_acb_poly_overlaps(y, len, yayb, len))
        {
            flint_printf("FAIL\n\n");
            flint_printf("len = %wd, k = %wu, n = %wu, m1 = %wu, m2 = %wu, m3 = %wu\n\n", len, k, n, m1, m2, m3);
            flint_printf("x = "); acb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); _acb_vec_printn(y, len, 100, 0); flint_printf("\n\n");
            flint_printf("ya = "); _acb_vec_printn(ya, len, 100, 0); flint_printf("\n\n");
            flint_printf("yb = "); _acb_vec_printn(yb, len, 100, 0); flint_printf("\n\n");
            flint_printf("yayb = "); _acb_vec_printn(yayb, len, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(xk);
        _acb_vec_clear(y, len);
        _acb_vec_clear(ya, len);
        _acb_vec_clear(yb, len);
        _acb_vec_clear(yayb, len);
    }

    TEST_FUNCTION_END(state);
}
#undef rising_algorithm
