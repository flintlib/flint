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
#include "arb_poly.h"

TEST_FUNCTION_START(arb_rising2_ui, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, u, v, u2, v2;
        fmpz *f;
        arb_ptr g;
        ulong n;
        slong i, prec;

        arb_init(a);
        arb_init(u);
        arb_init(v);
        arb_init(u2);
        arb_init(v2);

        arb_randtest(a, state, 1 + n_randint(state, 4000), 10);
        arb_randtest(u, state, 1 + n_randint(state, 4000), 10);
        arb_randtest(v, state, 1 + n_randint(state, 4000), 10);
        n = n_randint(state, 120);

        f = _fmpz_vec_init(n + 1);
        g = _arb_vec_init(n + 1);

        prec = 2 + n_randint(state, 4000);
        arb_rising2_ui(u, v, a, n, prec);

        arith_stirling_number_1u_vec(f, n, n + 1);
        for (i = 0; i <= n; i++)
            arb_set_fmpz(g + i, f + i);
        _arb_poly_evaluate(u2, g, n + 1, a, prec);

        _arb_poly_derivative(g, g, n + 1, prec);
        _arb_poly_evaluate(v2, g, n, a, prec);

        if (!arb_overlaps(u, u2) || !arb_overlaps(v, v2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n", n);
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("u = "); arb_printd(u, 15); flint_printf("\n\n");
            flint_printf("u2 = "); arb_printd(u2, 15); flint_printf("\n\n");
            flint_printf("v = "); arb_printd(v, 15); flint_printf("\n\n");
            flint_printf("v2 = "); arb_printd(v2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_set(u2, a);
        arb_rising2_ui(u2, v, u2, n, prec);

        if (!arb_equal(u2, u))
        {
            flint_printf("FAIL: aliasing 1\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("u = "); arb_printd(u, 15); flint_printf("\n\n");
            flint_printf("u2 = "); arb_printd(u2, 15); flint_printf("\n\n");
            flint_printf("n = %wu\n", n);
            flint_abort();
        }

        arb_set(v2, a);
        arb_rising2_ui(u, v2, v2, n, prec);

        if (!arb_equal(v2, v))
        {
            flint_printf("FAIL: aliasing 2\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("v = "); arb_printd(v, 15); flint_printf("\n\n");
            flint_printf("v2 = "); arb_printd(v2, 15); flint_printf("\n\n");
            flint_printf("n = %wu\n", n);
            flint_abort();
        }

        arb_clear(a);
        arb_clear(u);
        arb_clear(v);
        arb_clear(u2);
        arb_clear(v2);
        _fmpz_vec_clear(f, n + 1);
        _arb_vec_clear(g, n + 1);
    }

    TEST_FUNCTION_END(state);
}
