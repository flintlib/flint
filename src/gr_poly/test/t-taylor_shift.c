/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("taylor_shift....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t F, Fc, Fcd1, Fcd2;
        gr_ptr c, d, cd;

        /* check F(x+c)(x+d) = F(x+c+d) */
        gr_ctx_init_random(ctx, state);

        gr_poly_init(F, ctx);
        gr_poly_init(Fc, ctx);
        gr_poly_init(Fcd1, ctx);
        gr_poly_init(Fcd2, ctx);
        GR_TMP_INIT3(c, d, cd, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(F, state, 1 + n_randint(state, 10), ctx);
        status |= gr_poly_randtest(Fc, state, 1 + n_randint(state, 10), ctx);
        status |= gr_randtest(c, state, ctx);
        status |= gr_randtest(d, state, ctx);
        status |= gr_add(cd, c, d, ctx);

        status |= gr_poly_taylor_shift(Fc, F, c, ctx);
        status |= gr_poly_taylor_shift(Fcd1, Fc, d, ctx);
        status |= gr_poly_set(Fcd2, F, ctx); /* also test aliasing */
        status |= gr_poly_taylor_shift(Fcd2, Fcd2, cd, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(Fcd1, Fcd2, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("c = "); gr_println(c, ctx);
            flint_printf("d = "); gr_println(d, ctx);
            flint_printf("cd = "); gr_println(cd, ctx);
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("Fc = "); gr_poly_print(Fc, ctx); flint_printf("\n");
            flint_printf("Fcd1 = "); gr_poly_print(Fcd1, ctx); flint_printf("\n");
            flint_printf("Fcd2 = "); gr_poly_print(Fcd2, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(F, ctx);
        gr_poly_clear(Fc, ctx);
        gr_poly_clear(Fcd1, ctx);
        gr_poly_clear(Fcd2, ctx);
        GR_TMP_CLEAR3(c, d, cd, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
