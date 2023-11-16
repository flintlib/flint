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

TEST_FUNCTION_START(ca_conj, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, z;
        truth_t equal;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(z, ctx);

        ca_randtest_special(x, state, 5, 5, ctx);

        switch (n_randint(state, 3))
        {
            case 0: ca_conj(y, x, ctx); break;
            case 1: ca_conj_deep(y, x, ctx); break;
            case 2: ca_conj_shallow(y, x, ctx); break;
        }

        switch (n_randint(state, 3))
        {
            case 0: ca_conj(z, y, ctx); break;
            case 1: ca_conj_deep(z, y, ctx); break;
            case 2: ca_conj_shallow(z, y, ctx); break;
        }

        equal = ca_check_equal(x, z, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(z, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
