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
test_atan_series(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong n;
    gr_poly_t a, b, c, d;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(a, ctx);
    gr_poly_init(b, ctx);
    gr_poly_init(c, ctx);
    gr_poly_init(d, ctx);

    n = n_randint(state, 5);

    GR_MUST_SUCCEED(gr_poly_randtest(a, state, 8, ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(b, state, 8, ctx));

    switch (which)
    {
        case 0:
            status |= gr_poly_atan_series(b, a, n, ctx);
            break;
        case 1:
            status |= gr_poly_set(b, a, ctx);
            status |= gr_poly_atan_series(b, b, n, ctx);
            break;
        case 2:
            status |= gr_poly_atanh_series(b, a, n, ctx);
            break;
        case 3:
            status |= gr_poly_set(b, a, ctx);
            status |= gr_poly_atanh_series(b, b, n, ctx);
            break;
        default:
            flint_abort();
    }

    if (status == GR_SUCCESS)
    {
        /* Check 2 atan(a) = atan(2a/(1-a^2)) + C */
        /*       2 atanh(a) = atanh(2a/(1+a^2)) + C */
        status |= gr_poly_mullow(c, a, a, n, ctx);
        status |= gr_poly_one(d, ctx);
        if (which == 0 || which == 1)
            status |= gr_poly_sub(c, d, c, ctx);
        else
            status |= gr_poly_add(c, d, c, ctx);

        status |= gr_poly_add(d, a, a, ctx);  /* todo: should have mul_2exp / mul_two */
        status |= gr_poly_div_series(c, d, c, n, ctx);

        if (which == 0 || which == 1)
            status |= gr_poly_atan_series(c, c, n, ctx);
        else
            status |= gr_poly_atanh_series(c, c, n, ctx);

        status |= gr_poly_add(d, b, b, ctx);

        /* todo: also check the first coefficient */
        status |= gr_poly_set_coeff_si(c, 0, 0, ctx);
        status |= gr_poly_set_coeff_si(d, 0, 0, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(c, d, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("which = %d, n = %wd\n", which, n);
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n\n");
            flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n\n");
            flint_printf("c = "); gr_poly_print(c, ctx); flint_printf("\n\n");
            flint_printf("d = "); gr_poly_print(d, ctx); flint_printf("\n\n");
            flint_abort();
        }
    }

    gr_poly_clear(a, ctx);
    gr_poly_clear(b, ctx);
    gr_poly_clear(c, ctx);
    gr_poly_clear(d, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_atan_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_atan_series(state, n_randint(state, 4));
    }

    TEST_FUNCTION_END(state);
}
