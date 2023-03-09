/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_poly.h"

extern gr_static_method_table _ca_methods;

int main()
{
    slong iter;
    slong count_success = 0, count_unable = 0, count_domain = 0;
    flint_rand_t state;

    flint_printf("adjugate....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        int status = GR_SUCCESS;
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A, B, C;
        gr_ptr d, e;

        gr_ctx_init_random(ctx, state);

        /* Hack: avoid because slow */
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        n = n_randint(state, 5);
        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);

        GR_TMP_INIT2(d, e, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        status |= gr_mat_adjugate_cofactor(B, d, A, ctx);
        status |= gr_mat_adjugate_charpoly(C, e, A, ctx);

        if (status == GR_SUCCESS && (gr_mat_equal(B, C, ctx) == T_FALSE || gr_equal(d, e, ctx) == T_FALSE))
        {
            flint_printf("FAIL\n");
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "), gr_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "), gr_mat_print(C, ctx); flint_printf("\n");
            flint_printf("d = "), gr_print(d, ctx); flint_printf("\n");
            flint_printf("e = "), gr_print(e, ctx); flint_printf("\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);

        GR_TMP_CLEAR2(d, e, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" [%wd success, %wd domain, %wd unable] PASS\n", count_success, count_domain, count_unable);
    return 0;
}
