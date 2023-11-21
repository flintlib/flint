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

TEST_FUNCTION_START(arb_chebyshev_u_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d, e;
        ulong n;
        slong prec;

        n = n_randtest(state);
        prec = 2 + n_randint(state, 300);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);

        arb_randtest_precise(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest_precise(c, state, 1 + n_randint(state, 300), 5);

        arb_sin_cos(d, b, a, prec);
        arb_chebyshev_u_ui(c, n, b, prec);
        arb_mul(d, c, d, prec);

        if (n == LIMB_ONES)
            arb_mul_2exp_si(e, a, FLINT_BITS);
        else
            arb_mul_ui(e, a, n + 1, prec);

        arb_sin(e, e, prec);

        if (!arb_overlaps(d, e))
        {
            flint_printf("FAIL: sin(a)*U_n(cos(a)) = sin((n+1)a)\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 15); flint_printf("\n\n");
            flint_printf("e = "); arb_printd(e, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_chebyshev_u_ui(b, n, b, prec);

        if (!arb_equal(b, c))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest(b, state, 1 + n_randint(state, 300), 5);
        arb_randtest(c, state, 1 + n_randint(state, 300), 5);

        arb_chebyshev_u2_ui(b, c, n, a, prec);
        arb_chebyshev_u_ui(d, n, a, prec);

        if (!arb_overlaps(b, d))
        {
            flint_printf("FAIL: U_n\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        if (n == 0)
            arb_zero(d);
        else
            arb_chebyshev_u_ui(d, n - 1, a, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: U_{n-1}\n\n");
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
        arb_clear(e);
    }

    TEST_FUNCTION_END(state);
}
