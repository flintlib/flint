/*
    Copyright (C) 2026 Maria Neagoie

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_ore_poly.h"

TEST_FUNCTION_START(gr_ore_poly_divrem, state)
{
    slong iter;
    // slong success = 0, domain = 0, unable = 0, combined = 0;
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_ore_poly_ctx_t ore_ctx;
        gr_ore_poly_t A, B, Q, R, QBR;

        gr_ore_poly_ctx_init_randtest2(ctx, ore_ctx, state);

        gr_ore_poly_init(A, ore_ctx);
        gr_ore_poly_init(B, ore_ctx);
        gr_ore_poly_init(Q, ore_ctx);
        gr_ore_poly_init(R, ore_ctx);
        gr_ore_poly_init(QBR, ore_ctx);

        status |= gr_ore_poly_randtest(A, state, 1 + n_randint(state, 6), ore_ctx);
        status |= gr_ore_poly_randtest(B, state, 1 + n_randint(state, 6), ore_ctx);
        status |= gr_ore_poly_randtest(Q, state, 1 + n_randint(state, 6), ore_ctx);
        status |= gr_ore_poly_randtest(R, state, 1 + n_randint(state, 6), ore_ctx);

        if (n_randint(state, 3) == 0)
        {
            status |= gr_ore_poly_mul(A, A, B, ore_ctx);
            status |= gr_ore_poly_add(A, A, R, ore_ctx);
        }

        if (status == GR_SUCCESS)
        {
            /* test aliasing */
            switch (n_randint(state, 5))
            {
                case 0:
                    status |= gr_ore_poly_set(Q, A, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, Q, B, ore_ctx);
                    break;
                case 1:
                    status |= gr_ore_poly_set(R, A, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, R, B, ore_ctx);
                    break;
                case 2:
                    status |= gr_ore_poly_set(Q, B, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, A, Q, ore_ctx);
                    break;
                case 3:
                    status |= gr_ore_poly_set(R, B, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, A, R, ore_ctx);
                    break;
                default:
                    status |= gr_ore_poly_divrem(Q, R, A, B, ore_ctx);
                    break;
            }

            if (status == GR_SUCCESS)
            {
                status |= gr_ore_poly_mul(QBR, Q, B, ore_ctx);
                status |= gr_ore_poly_add(QBR, QBR, R, ore_ctx);

                if (status == GR_SUCCESS && gr_ore_poly_equal(QBR, A, ore_ctx) == T_FALSE)
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("A = "); gr_ore_poly_print(A, ore_ctx); flint_printf("\n");
                    flint_printf("B = "); gr_ore_poly_print(B, ore_ctx); flint_printf("\n");
                    flint_printf("Q = "); gr_ore_poly_print(Q, ore_ctx); flint_printf("\n");
                    flint_printf("R = "); gr_ore_poly_print(R, ore_ctx); flint_printf("\n");
                    flint_printf("Q*B + R = "); gr_ore_poly_print(QBR, ore_ctx); flint_printf("\n");
                    flint_abort();
                }
                // success++;
            }
        }

        /*if (status == GR_DOMAIN)
            domain++;

        if (status == GR_UNABLE)
            unable++;

        if (status == 3)
            combined++;*/

        gr_ore_poly_clear(A, ore_ctx);
        gr_ore_poly_clear(B, ore_ctx);
        gr_ore_poly_clear(Q, ore_ctx);
        gr_ore_poly_clear(R, ore_ctx);
        gr_ore_poly_clear(QBR, ore_ctx);

        gr_ore_poly_ctx_clear(ore_ctx);
        gr_ctx_clear(ctx);
    }
    /*flint_printf("GR_SUCCESS = %d\n", success);
    flint_printf("GR_DOMAIN = %d\n", domain);
    flint_printf("GR_UNABLE = %d\n", unable);
    flint_printf("Combined (3) = %d\n", combined);*/
    TEST_FUNCTION_END(state);
}
