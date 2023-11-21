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

/* Defined in t-pow_series_fmpq.c, t-pow_series_ui.c and t-pow_ui.c */
#define test test_pow_series_ui
int
test(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    /* Test A^n * A^m = A^(n+m) */
    gr_poly_t A, An, Am, AnAm, Anm;
    ulong n, m;
    slong len;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(A, ctx);
    gr_poly_init(An, ctx);
    gr_poly_init(Am, ctx);
    gr_poly_init(AnAm, ctx);
    gr_poly_init(Anm, ctx);

    n = n_randint(state, 5);
    m = n_randint(state, 5);
    len = n_randint(state, 10);

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 8, ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(An, state, 8, ctx));

    if (which == 0)
    {
        status |= gr_poly_pow_series_ui(An, A, n, len + n_randint(state, 5), ctx);
        status |= gr_poly_set(Am, A, ctx);
        status |= gr_poly_pow_series_ui(Am, Am, m, len + n_randint(state, 5), ctx);
        status |= gr_poly_mullow(AnAm, An, Am, len, ctx);
        status |= gr_poly_pow_series_ui(Anm, A, n + m, len, ctx);
    }
    else
    {
        status |= gr_poly_pow_series_ui_binexp(An, A, n, len + n_randint(state, 5), ctx);
        status |= gr_poly_set(Am, A, ctx);
        status |= gr_poly_pow_series_ui_binexp(Am, Am, m, len + n_randint(state, 5), ctx);
        status |= gr_poly_mullow(AnAm, An, Am, len, ctx);
        status |= gr_poly_pow_series_ui_binexp(Anm, A, n + m, len, ctx);
    }

    if (status == GR_SUCCESS && gr_poly_equal(AnAm, Anm, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n\n");
        gr_ctx_println(ctx);
        flint_printf("n = %wu, m = %wu, len = %wd\n", n, m, len);
        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
        flint_printf("An = "); gr_poly_print(An, ctx); flint_printf("\n");
        flint_printf("Am = "); gr_poly_print(Am, ctx); flint_printf("\n");
        flint_printf("AnAm = "); gr_poly_print(AnAm, ctx); flint_printf("\n");
        flint_printf("Anm = "); gr_poly_print(Anm, ctx); flint_printf("\n");
        flint_abort();
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(An, ctx);
    gr_poly_clear(Am, ctx);
    gr_poly_clear(AnAm, ctx);
    gr_poly_clear(Anm, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_pow_series_ui, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test(state, n_randint(state, 2));
    }

    TEST_FUNCTION_END(state);
}
#undef test
