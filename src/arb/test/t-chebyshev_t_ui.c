/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_chebyshev_t_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        ulong n;
        slong prec;

        n = n_randtest(state);
        prec = 2 + n_randint(state, 300);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest(c, state, 1 + n_randint(state, 300), 5);

        arb_cos(b, a, prec);
        arb_chebyshev_t_ui(c, n, b, prec);

        arb_mul_ui(d, a, n, prec);
        arb_cos(d, d, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: c = T_n(cos(a)) = d = cos(n*a)\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("d = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_chebyshev_t_ui(b, n, b, prec);

        if (!arb_equal(b, c))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest(b, state, 1 + n_randint(state, 300), 5);
        arb_randtest(c, state, 1 + n_randint(state, 300), 5);

        arb_chebyshev_t2_ui(b, c, n, a, prec);
        arb_chebyshev_t_ui(d, n, a, prec);

        if (!arb_overlaps(b, d))
        {
            flint_printf("FAIL: T_n\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        if (n == 0)
            arb_set(d, a);
        else
            arb_chebyshev_t_ui(d, n - 1, a, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: T_{n-1}\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
