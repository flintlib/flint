/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"

TEST_FUNCTION_START(ca_gamma, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, gx, gy, gxgy;
        acb_t gax, gay, agx, agy, agxgy, gaxgay;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(gx, ctx);
        ca_init(gy, ctx);
        ca_init(gxgy, ctx);
        acb_init(gax);
        acb_init(gay);
        acb_init(agx);
        acb_init(agy);
        acb_init(agxgy);
        acb_init(gaxgay);

        ca_randtest(x, state, 3, 5, ctx);
        if (n_randint(state, 2))
            ca_randtest(y, state, 3, 5, ctx);
        else
            ca_randtest_rational(y, state, 5, ctx);
        if (n_randint(state, 2))
            ca_add(y, y, x, ctx);

        ca_gamma(gx, x, ctx);
        ca_gamma(gy, y, ctx);
        ca_add(gxgy, gx, gy, ctx);

        ca_get_acb(agx, gx, 53, ctx);
        ca_get_acb(agy, gy, 53, ctx);
        ca_get_acb(agxgy, gxgy, 53, ctx);

        ca_get_acb(gax, x, 53, ctx);
        acb_gamma(gax, gax, 53);
        ca_get_acb(gay, y, 53, ctx);
        acb_gamma(gay, gay, 53);
        acb_add(gaxgay, gax, gay, 53);

        if (!acb_overlaps(agx, gax) || !acb_overlaps(agy, gay) ||
            !acb_overlaps(agxgy, gaxgay))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); ca_print(x, ctx); printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); printf("\n\n");
            flint_printf("gx = "); ca_print(gx, ctx); printf("\n\n");
            flint_printf("gy = "); ca_print(gy, ctx); printf("\n\n");
            flint_printf("gxgy = "); ca_print(gxgy, ctx); printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(gx, ctx);
        ca_clear(gy, ctx);
        ca_clear(gxgy, ctx);
        acb_clear(gax);
        acb_clear(gay);
        acb_clear(agx);
        acb_clear(agy);
        acb_clear(agxgy);
        acb_clear(gaxgay);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
