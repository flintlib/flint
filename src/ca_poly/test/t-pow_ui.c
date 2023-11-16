/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_poly.h"

TEST_FUNCTION_START(ca_poly_pow_ui, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, An, Am, AnAm, Anm;
        ulong n, m;

        /* Test A^n * A^m = A^(n+m) */
        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(An, ctx);
        ca_poly_init(Am, ctx);
        ca_poly_init(AnAm, ctx);
        ca_poly_init(Anm, ctx);

        ca_poly_randtest(A, state, 5, 1, 5, ctx);
        n = n_randint(state, 5);
        m = n_randint(state, 5);

        ca_poly_pow_ui(An, A, n, ctx);
        ca_poly_set(Am, A, ctx);
        ca_poly_pow_ui(Am, Am, m, ctx);
        ca_poly_mul(AnAm, An, Am, ctx);
        ca_poly_pow_ui(Anm, A, n + m, ctx);

        if (ca_poly_check_equal(AnAm, Anm, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
            flint_printf("An = "); ca_poly_print(An, ctx); flint_printf("\n");
            flint_printf("Am = "); ca_poly_print(Am, ctx); flint_printf("\n");
            flint_printf("AnAm = "); ca_poly_print(AnAm, ctx); flint_printf("\n");
            flint_printf("Anm = "); ca_poly_print(Anm, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_poly_clear(A, ctx);
        ca_poly_clear(An, ctx);
        ca_poly_clear(Am, ctx);
        ca_poly_clear(AnAm, ctx);
        ca_poly_clear(Anm, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
