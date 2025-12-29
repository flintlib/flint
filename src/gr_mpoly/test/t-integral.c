/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_integral, state)
{
    slong iter;

    for (iter = 0; iter < 10000; iter++)
    {
        int status;
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t F, G, Gprime;
        slong len;
        flint_bitcnt_t exp_bits;
        slong var;

        /* check derivative(integral(f)) == f */

        gr_ctx_init_random(cctx, state);
        gr_mpoly_ctx_init_rand(ctx, state, cctx, 10);
        /* ensure we have a variable */
        while (GR_MPOLY_NVARS(ctx) == 0)
        {
            gr_mpoly_ctx_clear(ctx);
            gr_mpoly_ctx_init_rand(ctx, state, cctx, 10);
        }

        var = n_randint(state, GR_MPOLY_NVARS(ctx));

        gr_mpoly_init(F, ctx);
        gr_mpoly_init(G, ctx);
        gr_mpoly_init(Gprime, ctx);

        status = GR_SUCCESS;

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;

        status |= gr_mpoly_randtest_bits(F, state, len, exp_bits, ctx);

        if (n_randint(state, 2))
        {
            status |= gr_mpoly_integral(G, F, var, ctx);
        }
        else
        {
            status |= gr_mpoly_set(G, F, ctx);
            status |= gr_mpoly_integral(G, G, var, ctx);
        }

        if (n_randint(state, 2))
        {
            status |= gr_mpoly_derivative(Gprime, G, var, ctx);
        }
        else
        {
            status |= gr_mpoly_set(Gprime, G, ctx);
            status |= gr_mpoly_derivative(Gprime, Gprime, var, ctx);
        }

        if (status == GR_SUCCESS && gr_mpoly_equal(F, Gprime, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("var = %ld\n", var);
            flint_printf("F = "); gr_mpoly_print_pretty(F, ctx); flint_printf("\n");
            flint_printf("G = "); gr_mpoly_print_pretty(G, ctx); flint_printf("\n");
            flint_printf("Gprime = "); gr_mpoly_print_pretty(Gprime, ctx); flint_printf("\n");
            status |= gr_mpoly_sub(Gprime, Gprime, F, ctx);
            flint_printf("Gprime - F = "); gr_mpoly_print_pretty(Gprime, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_mpoly_clear(F, ctx);
        gr_mpoly_clear(G, ctx);
        gr_mpoly_clear(Gprime, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}

