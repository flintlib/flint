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

    flint_printf("acos...");
    fflush(stdout);

    flint_randinit(state);

    {
        ca_ctx_t ctx;
        ca_ctx_init(ctx);

        for (iter = 0; iter < 200 * calcium_test_multiplier(); iter++)
        {
            ca_t x, t1, t2, t3;

            if (n_randint(state, 10) == 0)
            {
                ca_ctx_clear(ctx);
                ca_ctx_init(ctx);
            }

            ca_init(x, ctx);
            ca_init(t1, ctx);
            ca_init(t2, ctx);
            ca_init(t3, ctx);
     
            ca_randtest_special(x, state, 5, 5, ctx);

            ca_acos_direct(t1, x, ctx);
            ca_acos_logarithm(t2, x, ctx);

            if (ca_check_equal(t1, t2, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                flint_printf("t1 = "); ca_print(t1, ctx); flint_printf("\n\n");
                flint_printf("t2 = "); ca_print(t2, ctx); flint_printf("\n\n");
                flint_abort();
            }

            if (n_randint(state, 2))
                ca_cos(t3, t2, ctx);
            else
                ca_cos(t3, t1, ctx);

            if (ca_check_is_infinity(t3, ctx) != T_TRUE && ca_check_is_undefined(t3, ctx) != T_TRUE)
            {
                if (ca_check_equal(t3, x, ctx) == T_FALSE)
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
                    flint_printf("t2 = "); ca_print(t2, ctx); flint_printf("\n\n");
                    flint_printf("t3 = "); ca_print(t3, ctx); flint_printf("\n\n");
                    flint_abort();
                }
            }

            ca_clear(x, ctx);
            ca_clear(t1, ctx);
            ca_clear(t2, ctx);
            ca_clear(t3, ctx);
        }

        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
