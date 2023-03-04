/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("evaluate....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t F, G, FG;
        gr_ptr x, Fx, Gx, FxGx, FGx;

        /* Test F(x) + G(x) = (F + G)(x) */
        gr_ctx_init_random(ctx, state);

        gr_poly_init(F, ctx);
        gr_poly_init(G, ctx);
        gr_poly_init(FG, ctx);
        GR_TMP_INIT5(x, Fx, Gx, FxGx, FGx, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(F, state, 1 + n_randint(state, 10), ctx);
        status |= gr_poly_randtest(G, state, 1 + n_randint(state, 10), ctx);
        status |= gr_randtest(x, state, ctx);
        status |= gr_randtest(Fx, state, ctx);

        status |= gr_poly_add(FG, F, G, ctx);

        status |= gr_poly_evaluate(Fx, F, x, ctx);
        status |= gr_set(Gx, x, ctx); /* test aliasing */
        status |= gr_poly_evaluate(Gx, G, Gx, ctx);
        status |= gr_add(FxGx, Fx, Gx, ctx);
        status |= gr_poly_evaluate(FGx, FG, x, ctx);

        if (status == GR_SUCCESS && gr_equal(FxGx, FGx, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
            flint_printf("x = "); gr_print(x, ctx); flint_printf("\n");
            flint_printf("FxGx = "); gr_print(FxGx, ctx); flint_printf("\n");
            flint_printf("FGx = "); gr_print(FGx, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(F, ctx);
        gr_poly_clear(G, ctx);
        gr_poly_clear(FG, ctx);
        GR_TMP_CLEAR5(x, Fx, Gx, FxGx, FGx, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
