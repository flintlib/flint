/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_rising2_ui, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, u, v, u2, v2;
        fmpz *f;
        acb_ptr g;
        ulong n;
        slong i, prec;

        acb_init(a);
        acb_init(u);
        acb_init(v);
        acb_init(u2);
        acb_init(v2);

        acb_randtest(a, state, 1 + n_randint(state, 4000), 10);
        acb_randtest(u, state, 1 + n_randint(state, 4000), 10);
        acb_randtest(v, state, 1 + n_randint(state, 4000), 10);
        n = n_randint(state, 120);

        f = _fmpz_vec_init(n + 1);
        g = _acb_vec_init(n + 1);

        prec = 2 + n_randint(state, 4000);
        acb_rising2_ui(u, v, a, n, prec);

        arith_stirling_number_1u_vec(f, n, n + 1);
        for (i = 0; i <= n; i++)
            acb_set_fmpz(g + i, f + i);
        _acb_poly_evaluate(u2, g, n + 1, a, prec);

        _acb_poly_derivative(g, g, n + 1, prec);
        _acb_poly_evaluate(v2, g, n, a, prec);

        if (!acb_overlaps(u, u2) || !acb_overlaps(v, v2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n", n);
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("u = "); acb_printd(u, 15); flint_printf("\n\n");
            flint_printf("u2 = "); acb_printd(u2, 15); flint_printf("\n\n");
            flint_printf("v = "); acb_printd(v, 15); flint_printf("\n\n");
            flint_printf("v2 = "); acb_printd(v2, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_set(u2, a);
        acb_rising2_ui(u2, v, u2, n, prec);

        if (!acb_equal(u2, u))
        {
            flint_printf("FAIL: aliasing 1\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("u = "); acb_printd(u, 15); flint_printf("\n\n");
            flint_printf("u2 = "); acb_printd(u2, 15); flint_printf("\n\n");
            flint_printf("n = %wu\n", n);
            flint_abort();
        }

        acb_set(v2, a);
        acb_rising2_ui(u, v2, v2, n, prec);

        if (!acb_equal(v2, v))
        {
            flint_printf("FAIL: aliasing 2\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("v = "); acb_printd(v, 15); flint_printf("\n\n");
            flint_printf("v2 = "); acb_printd(v2, 15); flint_printf("\n\n");
            flint_printf("n = %wu\n", n);
            flint_abort();
        }

        acb_clear(a);
        acb_clear(u);
        acb_clear(v);
        acb_clear(u2);
        acb_clear(v2);
        _fmpz_vec_clear(f, n + 1);
        _acb_vec_clear(g, n + 1);
    }

    TEST_FUNCTION_END(state);
}
