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

TEST_FUNCTION_START(ca_poly_evaluate_horner, state)
{
    slong iter;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t F, G, FG;
        ca_t x, Fx, Gx, FxGx, FGx;

        /* Test F(x) + G(x) = (F + G)(x) */
        ca_ctx_init(ctx);

        ca_poly_init(F, ctx);
        ca_poly_init(G, ctx);
        ca_poly_init(FG, ctx);
        ca_init(x, ctx);
        ca_init(Fx, ctx);
        ca_init(Gx, ctx);
        ca_init(FxGx, ctx);
        ca_init(FGx, ctx);

        ca_poly_randtest(F, state, 8, 1, 5, ctx);
        ca_poly_randtest(G, state, 8, 1, 5, ctx);
        ca_randtest(x, state, 1, 5, ctx);
        ca_randtest(Fx, state, 1, 5, ctx);

        ca_poly_add(FG, F, G, ctx);

        ca_poly_evaluate_horner(Fx, F, x, ctx);
        ca_set(Gx, x, ctx); /* test aliasing */
        ca_poly_evaluate_horner(Gx, G, Gx, ctx);
        ca_add(FxGx, Fx, Gx, ctx);
        ca_poly_evaluate_horner(FGx, FG, x, ctx);

        if (ca_check_equal(FxGx, FGx, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("F = "); ca_poly_print(F, ctx); flint_printf("\n");
            flint_printf("G = "); ca_poly_print(G, ctx); flint_printf("\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n");
            flint_printf("FxGx = "); ca_print(FxGx, ctx); flint_printf("\n");
            flint_printf("FGx = "); ca_print(FGx, ctx); flint_printf("\n");
            flint_abort();
        }

        ca_poly_clear(F, ctx);
        ca_poly_clear(G, ctx);
        ca_poly_clear(FG, ctx);
        ca_clear(x, ctx);
        ca_clear(Fx, ctx);
        ca_clear(Gx, ctx);
        ca_clear(FxGx, ctx);
        ca_clear(FGx, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
