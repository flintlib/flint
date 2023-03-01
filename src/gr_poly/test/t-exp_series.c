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

int
test_exp_series(flint_rand_t state)
{
    gr_ctx_t ctx;
    slong i, len1, len2, len3, n;
    gr_poly_t a, b, ab, fa, fb, fab, fafb;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(a, ctx);
    gr_poly_init(b, ctx);
    gr_poly_init(ab, ctx);
    gr_poly_init(fa, ctx);
    gr_poly_init(fb, ctx);
    gr_poly_init(fab, ctx);
    gr_poly_init(fafb, ctx);

    if (ctx->methods == _ca_methods)
    {
        len1 = n_randint(state, 5);
        len2 = n_randint(state, 5);
        len3 = n_randint(state, 5);
    }
    else
    {
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        len3 = n_randint(state, 20);
    }

    GR_MUST_SUCCEED(gr_poly_randtest(a, state, 1 + n_randint(state, 20), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(b, state, 1 + n_randint(state, 20), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(fa, state, 1 + n_randint(state, 20), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(fb, state, 1 + n_randint(state, 20), ctx));

    if (n_randint(state, 4) != 0)
        status |= gr_poly_set_coeff_si(a, 0, 0, ctx);
    if (n_randint(state, 4) != 0)
        status |= gr_poly_set_coeff_si(b, 0, 0, ctx);

    status |= gr_poly_add(ab, a, b, ctx);

    for (i = 0; i < 3 && status == GR_SUCCESS; i++)
    {
        gr_ptr x, y;

        if (i == 0)
        {
            x = a;
            y = fa;
            n = len1;
        }
        else if (i == 1)
        {
            x = b;
            y = fb;
            n = len2;
        }
        else
        {
            x = ab;
            y = fab;
            n = len3;
        }

        switch (n_randint(state, 8))
        {
            case 0:
                status |= gr_poly_exp_series(y, x, n, ctx);
                break;
            case 1:
                status |= gr_poly_set(y, x, ctx);
                status |= gr_poly_exp_series(y, y, n, ctx);
                break;
            case 2:
                status |= gr_poly_exp_series_basecase(y, x, n, ctx);
                break;
            case 3:
                status |= gr_poly_set(y, x, ctx);
                status |= gr_poly_exp_series_basecase(y, y, n, ctx);
                break;
            case 4:
                status |= gr_poly_exp_series_basecase_mul(y, x, n, ctx);
                break;
            case 5:
                status |= gr_poly_set(y, x, ctx);
                status |= gr_poly_exp_series_basecase_mul(y, y, n, ctx);
                break;
            case 6:
                status |= gr_poly_exp_series_newton(y, x, n, n_randint(state, 20), ctx);
                break;
            case 7:
                status |= gr_poly_set(y, x, ctx);
                status |= gr_poly_exp_series_newton(y, y, n, n_randint(state, 20), ctx);
                break;
            default:
                abort();
        }
    }

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_mullow(fafb, fa, fb, FLINT_MIN(FLINT_MIN(len1, len2), len3), ctx);
        status |= gr_poly_truncate(fab, FLINT_MIN(FLINT_MIN(len1, len2), len3), ctx);

        if (status == GR_SUCCESS && gr_poly_equal(fafb, fab, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("len1 = %wd, len2 = %wd, len3 = %wd\n", len1, len2, len3);
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n\n");
            flint_printf("ab = "); gr_poly_print(ab, ctx); flint_printf("\n\n");
            flint_printf("fa = "); gr_poly_print(fa, ctx); flint_printf("\n\n");
            flint_printf("fb = "); gr_poly_print(fb, ctx); flint_printf("\n\n");
            flint_printf("fab = "); gr_poly_print(fab, ctx); flint_printf("\n\n");
            flint_printf("fafb = "); gr_poly_print(fafb, ctx); flint_printf("\n\n");
            flint_abort();
        }
    }

    gr_poly_clear(a, ctx);
    gr_poly_clear(b, ctx);
    gr_poly_clear(ab, ctx);
    gr_poly_clear(fa, ctx);
    gr_poly_clear(fb, ctx);
    gr_poly_clear(fab, ctx);
    gr_poly_clear(fafb, ctx);

    gr_ctx_clear(ctx);

    return status;
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("exp_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        test_exp_series(state);
    }

    flint_randclear(state);
    flint_cleanup_master();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
