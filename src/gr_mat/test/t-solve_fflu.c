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

    flint_printf("solve_fflu....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, X, B, AX;
        slong n, c;
        int status = GR_SUCCESS;

        n = n_randint(state, 6);
        c = n_randint(state, 6);

        /* Hack: avoid because slow */
        gr_ctx_init_random(ctx, state);
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, c, ctx);
        gr_mat_init(X, n, c, ctx);
        gr_mat_init(AX, n, c, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_randtest(X, state, ctx);
        status |= gr_mat_randtest(B, state, ctx);

        status |= gr_mat_nonsingular_solve_fflu(X, A, B, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_mat_mul(AX, A, X, ctx);

            if (status == GR_SUCCESS && gr_mat_equal(AX, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n");
                gr_ctx_println(ctx);
                flint_printf("A = \n"); gr_mat_print(A, ctx); flint_printf("\n\n");
                flint_printf("B = \n"); gr_mat_print(B, ctx); flint_printf("\n\n");
                flint_printf("X = \n"); gr_mat_print(X, ctx); flint_printf("\n\n");
                flint_printf("AX = \n"); gr_mat_print(AX, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }
        else if ((status & GR_DOMAIN) && !(status & GR_UNABLE) && gr_ctx_is_integral_domain(ctx) == T_TRUE)
        {
            gr_ptr det;
            int status2;

            GR_TMP_INIT(det, ctx);
            status2 = gr_mat_det_berkowitz(det, A, ctx);

            if (status2 == GR_SUCCESS && gr_is_invertible(det, ctx) == T_TRUE)
            {
                flint_printf("FAIL (singular)\n\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("det = "); gr_print(det, ctx); flint_printf("\n");
                flint_abort();
            }

            GR_TMP_CLEAR(det, ctx);
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(X, ctx);
        gr_mat_clear(AX, ctx);
        gr_ctx_clear(ctx);

    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf(" [%wd success, %wd domain, %wd unable] PASS\n", count_success, count_domain, count_unable);
    return EXIT_SUCCESS;
}
