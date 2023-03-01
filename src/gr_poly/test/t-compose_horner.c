/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

extern gr_static_method_table _ca_methods;

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("compose_horner....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t F, G, FG;
        gr_poly_t x, Fx, Gx, FxGx, FGx;

        /* Test F(x) + G(x) = (F + G)(x) */
        gr_ctx_init_random(ctx, state);

        /* Hack: avoid because slow */
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        gr_poly_init(F, ctx);
        gr_poly_init(G, ctx);
        gr_poly_init(FG, ctx);
        gr_poly_init(x, ctx);
        gr_poly_init(Fx, ctx);
        gr_poly_init(Gx, ctx);
        gr_poly_init(FxGx, ctx);
        gr_poly_init(FGx, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(F, state, 1 + n_randint(state, 12), ctx);
        status |= gr_poly_randtest(G, state, 1 + n_randint(state, 12), ctx);
        status |= gr_poly_randtest(x, state, 1 + n_randint(state, 12), ctx);
        status |= gr_poly_randtest(Fx, state, 1 + n_randint(state, 12), ctx);

        status |= gr_poly_add(FG, F, G, ctx);

        status |= gr_poly_compose_horner(Fx, F, x, ctx);
        status |= gr_poly_set(Gx, x, ctx); /* test aliasing */
        status |= gr_poly_compose_horner(Gx, G, Gx, ctx);
        status |= gr_poly_add(FxGx, Fx, Gx, ctx);
        status |= gr_poly_compose_horner(FGx, FG, x, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(FxGx, FGx, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
            flint_printf("x = "); gr_poly_print(x, ctx); flint_printf("\n");
            flint_printf("Fx = "); gr_poly_print(Fx, ctx); flint_printf("\n");
            flint_printf("Gx = "); gr_poly_print(Gx, ctx); flint_printf("\n");
            flint_printf("FxGx = "); gr_poly_print(FxGx, ctx); flint_printf("\n");
            flint_printf("FGx = "); gr_poly_print(FGx, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(F, ctx);
        gr_poly_clear(G, ctx);
        gr_poly_clear(FG, ctx);
        gr_poly_clear(x, ctx);
        gr_poly_clear(Fx, ctx);
        gr_poly_clear(Gx, ctx);
        gr_poly_clear(FxGx, ctx);
        gr_poly_clear(FGx, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
