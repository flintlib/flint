/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

int
test_log_series(flint_rand_t state)
{
    gr_ctx_t ctx;
    slong n;
    gr_poly_t a, b, ab, fa, fb, fafb, fab;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(a, ctx);
    gr_poly_init(b, ctx);
    gr_poly_init(ab, ctx);
    gr_poly_init(fa, ctx);
    gr_poly_init(fb, ctx);
    gr_poly_init(fafb, ctx);
    gr_poly_init(fab, ctx);

    n = n_randint(state, 10);

    GR_MUST_SUCCEED(gr_poly_randtest(a, state, 10, ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(b, state, 10, ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(fa, state, 10, ctx));

    if (n_randint(state, 2))
        status |= gr_poly_set_coeff_ui(a, 0, 1, ctx);
    if (n_randint(state, 2))
        status |= gr_poly_set_coeff_ui(b, 0, 1, ctx);

    status |= gr_poly_log_series(fa, a, n, ctx);
    status |= gr_poly_log_series(fb, b, n, ctx);

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_add(fafb, fa, fb, ctx);
        status |= gr_poly_mullow(ab, a, b, n, ctx);
        status |= gr_poly_set(fab, ab, ctx);  /* also test aliasing */
        status |= gr_poly_log_series(fab, fab, n, ctx);

        /* todo: also test constant terms */
        status |= gr_poly_set_coeff_si(fab, 0, 0, ctx);
        status |= gr_poly_set_coeff_si(fafb, 0, 0, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(fab, fafb, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("n = %wd\n", n);
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n\n");
            flint_printf("fa = "); gr_poly_print(fa, ctx); flint_printf("\n\n");
            flint_printf("fb = "); gr_poly_print(fb, ctx); flint_printf("\n\n");
            flint_printf("ab = "); gr_poly_print(ab, ctx); flint_printf("\n\n");
            flint_printf("fafb = "); gr_poly_print(fafb, ctx); flint_printf("\n\n");
            flint_printf("fab = "); gr_poly_print(fab, ctx); flint_printf("\n\n");
            flint_abort();
        }

        status |= gr_poly_log_series(fb, a, n + n_randint(state, 3), ctx);
        status |= gr_poly_truncate(fb, fb, n, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(fa, fb, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("n = %wd\n", n);
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n\n");
            flint_printf("fa = "); gr_poly_print(fa, ctx); flint_printf("\n\n");
            flint_printf("fb = "); gr_poly_print(fb, ctx); flint_printf("\n\n");
            flint_abort();
        }
    }

    gr_poly_clear(a, ctx);
    gr_poly_clear(b, ctx);
    gr_poly_clear(ab, ctx);
    gr_poly_clear(fa, ctx);
    gr_poly_clear(fb, ctx);
    gr_poly_clear(fafb, ctx);
    gr_poly_clear(fab, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_log_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_log_series(state);
    }

    TEST_FUNCTION_END(state);
}
