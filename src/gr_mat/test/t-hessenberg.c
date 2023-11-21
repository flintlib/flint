/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_GR_FUNCTION_START(gr_mat_hessenberg, state, count_success, count_unable, count_domain)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status = GR_SUCCESS;
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A, B;
        gr_poly_t f, g;

        /* Hack: avoid because slow */
        gr_ctx_init_random(ctx, state);
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        n = n_randint(state, 7);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_hessenberg(B, A, ctx);

        status |= gr_mat_charpoly_berkowitz(f, A, ctx);

        if (n_randint(state, 2))
            status |= gr_mat_charpoly_from_hessenberg(g, B, ctx);
        else
            status |= gr_mat_charpoly_berkowitz(g, B, ctx);

        if (status == GR_SUCCESS &&
            (gr_poly_equal(f, g, ctx) == T_FALSE || gr_mat_is_hessenberg(B, ctx) == T_FALSE))
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_mat_print(B, ctx); flint_printf("\n");
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
