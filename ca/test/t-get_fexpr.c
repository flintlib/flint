/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_fexpr....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y;
        fexpr_t f;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        fexpr_init(f);

        ca_randtest_special(x, state, 5, 5, ctx);
        ca_get_fexpr(f, x, 0, ctx);

        if (!ca_set_fexpr(y, f, ctx))
        {
            flint_printf("FAIL: unable to parse\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("f = "); fexpr_print(f); flint_printf("\n\n");
            flint_abort();
        }

        if (ca_check_equal(x, y, ctx) == T_FALSE)
        {
            flint_printf("FAIL: not equal!\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("f = "); fexpr_print(f); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_abort();
        }

        if (ca_check_equal(x, y, ctx) == T_FALSE)
        {
            flint_printf("FAIL: not equal!\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("f = "); fexpr_print(f); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_abort();
        }

        if (ca_check_equal(x, y, ctx) != T_TRUE && !ca_is_unknown(x, ctx))
        {
            flint_printf("FAIL: not equal!\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("f = "); fexpr_print(f); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_abort();
        }

/*
        if (!ca_equal_repr(x, y, ctx))
        {
            flint_printf("Warning: not equal repr!\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("f = "); fexpr_print(f); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
        }
*/

        fexpr_clear(f);
        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
