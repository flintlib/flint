/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_chebyshev_u_ui, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d, e;
        ulong n;
        slong prec;

        n = n_randtest(state);
        prec = 2 + n_randint(state, 300);

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(e);

        acb_randtest_precise(a, state, 1 + n_randint(state, 300), 5);
        acb_randtest_precise(c, state, 1 + n_randint(state, 300), 5);

        acb_sin_cos(d, b, a, prec);
        acb_chebyshev_u_ui(c, n, b, prec);
        acb_mul(d, c, d, prec);

        if (n == LIMB_ONES)
            acb_mul_2exp_si(e, a, FLINT_BITS);
        else
            acb_mul_ui(e, a, n + 1, prec);

        acb_sin(e, e, prec);

        if (!acb_overlaps(d, e))
        {
            flint_printf("FAIL: sin(a)*U_n(cos(a)) = sin((n+1)a)\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 15); flint_printf("\n\n");
            flint_printf("e = "); acb_printd(e, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_chebyshev_u_ui(b, n, b, prec);

        if (!acb_equal(b, c))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_randtest(a, state, 1 + n_randint(state, 300), 5);
        acb_randtest(b, state, 1 + n_randint(state, 300), 5);
        acb_randtest(c, state, 1 + n_randint(state, 300), 5);

        acb_chebyshev_u2_ui(b, c, n, a, prec);
        acb_chebyshev_u_ui(d, n, a, prec);

        if (!acb_overlaps(b, d))
        {
            flint_printf("FAIL: U_n\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        if (n == 0)
            acb_zero(d);
        else
            acb_chebyshev_u_ui(d, n - 1, a, prec);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: U_{n-1}\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(e);
    }

    TEST_FUNCTION_END(state);
}
