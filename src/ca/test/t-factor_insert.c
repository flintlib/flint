/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"

TEST_FUNCTION_START(ca_factor_insert, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_factor_t fac;
        ca_t b, e, p, q, t;
        slong i, n;

        ca_ctx_init(ctx);
        ca_factor_init(fac, ctx);
        ca_init(b, ctx);
        ca_init(e, ctx);
        ca_init(p, ctx);
        ca_init(q, ctx);
        ca_init(t, ctx);

        /* insert repeated bases; exponents of equal bases must add up */
        ca_one(p, ctx);
        n = n_randint(state, 8);

        for (i = 0; i < n; i++)
        {
            switch (n_randint(state, 4))
            {
                case 0:
                    ca_set_ui(b, 2, ctx);
                    break;
                case 1:
                    ca_set_ui(b, 3, ctx);
                    break;
                case 2:
                    ca_pi(b, ctx);
                    break;
                default:
                    ca_sqrt_ui(b, 2, ctx);
                    break;
            }

            ca_set_ui(e, 1 + n_randint(state, 3), ctx);
            ca_factor_insert(fac, b, e, ctx);
            ca_pow(t, b, e, ctx);
            ca_mul(p, p, t, ctx);
        }

        ca_factor_get_ca(q, fac, ctx);

        if (ca_check_equal(p, q, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("p = "); ca_print(p, ctx); flint_printf("\n\n");
            flint_printf("q = "); ca_print(q, ctx); flint_printf("\n\n");
            flint_printf("fac = "); ca_factor_print(fac, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_factor_clear(fac, ctx);
        ca_clear(b, ctx);
        ca_clear(e, ctx);
        ca_clear(p, ctx);
        ca_clear(q, ctx);
        ca_clear(t, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
