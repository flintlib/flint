/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

void
rising_algorithm(arb_t res, const arb_t x, ulong n, ulong m, slong prec, int alg, int alias)
{
    if (alias)
    {
        arb_set(res, x);
        rising_algorithm(res, res, n, m, prec, alg, 0);
        return;
    }

    if (alg == 0)
        arb_hypgeom_rising_ui_rs(res, x, n, m, prec);
    else if (alg == 1)
        arb_hypgeom_rising_ui_forward(res, x, n, prec);
    else if (alg == 2)
        arb_hypgeom_rising_ui_bs(res, x, n, prec);
    else if (alg == 3)
        arb_hypgeom_rising_ui_rec(res, x, n, prec);
    else
        arb_hypgeom_rising_ui(res, x, n, prec);
}

TEST_FUNCTION_START(arb_hypgeom_rising_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, xk, y, ya, yb, yayb;
        ulong k, n, m1, m2, m3;
        slong prec;
        int alg1, alg2, alg3, alias1, alias2, alias3;

        prec = 2 + n_randint(state, 200);
        k = n_randint(state, 10);
        n = n_randint(state, 50);
        m1 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n + k, 1));
        m2 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(k, 1));
        m3 = n_randint(state, 2) ? 0 : 1 + n_randint(state, FLINT_MAX(n, 1));
        alg1 = n_randint(state, 5);
        alg2 = n_randint(state, 5);
        alg3 = n_randint(state, 5);
        alias1 = n_randint(state, 2);
        alias2 = n_randint(state, 2);
        alias3 = n_randint(state, 2);

        if (n_randint(state, 100) == 0)
            n += 100;

        arb_init(x);
        arb_init(xk);
        arb_init(y);
        arb_init(ya);
        arb_init(yb);
        arb_init(yayb);

        arb_randtest(x, state, prec, 10 + n_randint(state, 200));
        arb_add_ui(xk, x, k, prec);

        rising_algorithm(y, x, n + k, m1, prec, alg1, alias1);
        rising_algorithm(ya, x, k, m2, prec, alg2, alias2);
        rising_algorithm(yb, xk, n, m3, prec, alg3, alias3);
        arb_mul(yayb, ya, yb, prec);

        if (!arb_overlaps(y, yayb))
        {
            flint_printf("FAIL\n\n");
            flint_printf("k = %wu, n = %wu, m1 = %wu, m2 = %wu, m3 = %wu\n\n", k, n, m1, m2, m3);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("ya = "); arb_printn(ya, 100, 0); flint_printf("\n\n");
            flint_printf("yb = "); arb_printn(yb, 100, 0); flint_printf("\n\n");
            flint_printf("yayb = "); arb_printn(yayb, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(xk);
        arb_clear(y);
        arb_clear(ya);
        arb_clear(yb);
        arb_clear(yayb);
    }

    TEST_FUNCTION_END(state);
}
