/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_content_vars, state)
{
    slong i, j;

    /* check simple example */
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g;
        slong varl[2] = {0, 1};
        const char * vars[4] = {"x", "y", "z", "w"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, 11);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);

        nmod_mpoly_set_str_pretty(f, "x*y*(z+w)*(x*z+y*w+1+z*w)", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "z+w", vars, ctx);
        nmod_mpoly_content_vars(f, f, varl, 2, ctx);
        if (!nmod_mpoly_equal(f, g, ctx))
        {
            flint_printf("FAIL: check simple example\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, t;
        slong nvars, num_vars, len;
        ulong * exp_bounds;
        slong * vars;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        nvars = ctx->minfo->nvars;
        if (nvars < 1)
        {
            nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(t, ctx);

        exp_bounds = (ulong *) flint_malloc(nvars*sizeof(ulong));
        for (j = 0; j < nvars; j++)
            exp_bounds[j] = 1 + n_randint(state, 5);

        len = n_randint(state, 20);
        nmod_mpoly_randtest_bounds(f, state, len, exp_bounds, ctx);

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
        nmod_mpoly_randtest_bounds(t, state, len, exp_bounds, ctx);
        nmod_mpoly_mul(f, f, t, ctx);
        nmod_mpoly_repack_bits(f, f, f->bits + n_randint(state, FLINT_BITS), ctx);

        if (!nmod_mpoly_content_vars(g, f, vars, num_vars, ctx))
        {
            flint_printf("FAIL: check content could be computed\n");
            fflush(stdout);
            flint_abort();
        }

        if (nmod_mpoly_is_zero(g, ctx))
        {
            if (!nmod_mpoly_is_zero(f, ctx))
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
                if (nmod_mpoly_degree_si(g, vars[j], ctx) != 0)
                {
                    flint_printf("FAIL: content depends on a bad variable\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            if (!nmod_mpoly_divides(t, f, g, ctx))
            {
                flint_printf("FAIL: check content divides\n");
                fflush(stdout);
                flint_abort();
            }

            if (!nmod_mpoly_content_vars(t, t, vars, num_vars, ctx))
            {
                flint_printf("FAIL: check cofactor content could be computed\n");
                fflush(stdout);
                flint_abort();
            }

            if (!nmod_mpoly_is_one(t, ctx))
            {
                flint_printf("FAIL: check cofactor content is one\n");
                fflush(stdout);
                flint_abort();
            }
        }

        if (!nmod_mpoly_content_vars(f, f, vars, num_vars, ctx))
        {
            flint_printf("FAIL: check aliased content could be computed\n");
            fflush(stdout);
            flint_abort();
        }

        if (!nmod_mpoly_equal(f, g, ctx))
        {
            flint_printf("FAIL: check aliasing\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(exp_bounds);
        flint_free(vars);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
