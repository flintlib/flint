/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_content_vars, state)
{
    slong i, j;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, t;
        slong nvars, num_vars, len;
        ulong * exp_bounds;
        slong * vars;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        nvars = ctx->zctx->minfo->nvars;
        if (nvars < 1)
        {
            fmpq_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(t, ctx);

        exp_bounds = (ulong *) flint_malloc(nvars*sizeof(ulong));
        for (j = 0; j < nvars; j++)
            exp_bounds[j] = 1 + n_randint(state, 5);

        len = n_randint(state, 20);
        fmpq_mpoly_randtest_bounds(f, state, len, 200, exp_bounds, ctx);

        vars = (slong *) flint_malloc(nvars*sizeof(slong));
        for (j = 0; j < nvars; j++)
            vars[j] = j;

        for (j = 0; j < 2*nvars; j++)
        {
            slong k1 = n_randint(state, nvars);
            slong k2 = n_randint(state, nvars);
            FLINT_SWAP(slong, vars[k1], vars[k2]);
        }

        num_vars = 1 + n_randint(state, nvars);
        for (j = 0; j < num_vars; j++)
            exp_bounds[vars[j]] = 1;

        len = n_randint(state, 10);
        fmpq_mpoly_randtest_bounds(t, state, len, 200, exp_bounds, ctx);
        fmpq_mpoly_mul(f, f, t, ctx);
        fmpq_mpoly_repack_bits(f, f, f->zpoly->bits + n_randint(state, FLINT_BITS), ctx);

        if (!fmpq_mpoly_content_vars(g, f, vars, num_vars, ctx))
        {
            flint_printf("FAIL: check content could be computed\n");
            fflush(stdout);
            flint_abort();
        }

        if (fmpq_mpoly_is_zero(g, ctx))
        {
            if (!fmpq_mpoly_is_zero(f, ctx))
            {
                flint_printf("FAIL: check zero content\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            for (j = 0; j < num_vars; j++)
            {
                if (fmpq_mpoly_degree_si(g, vars[j], ctx) != 0)
                {
                    flint_printf("FAIL: content depends on a bad variable\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            if (!fmpq_mpoly_divides(t, f, g, ctx))
            {
                flint_printf("FAIL: check content divides\n");
                fflush(stdout);
                flint_abort();
            }

            if (!fmpq_mpoly_content_vars(t, t, vars, num_vars, ctx))
            {
                flint_printf("FAIL: check cofactor content could be computed\n");
                fflush(stdout);
                flint_abort();
            }

            if (!fmpq_mpoly_is_one(t, ctx))
            {
                flint_printf("FAIL: check cofactor content is one\n");
                fflush(stdout);
                flint_abort();
            }
        }

        if (!fmpq_mpoly_content_vars(f, f, vars, num_vars, ctx))
        {
            flint_printf("FAIL: check aliased content could be computed\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpq_mpoly_equal(f, g, ctx))
        {
            flint_printf("FAIL: check aliasing\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(exp_bounds);
        flint_free(vars);

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
