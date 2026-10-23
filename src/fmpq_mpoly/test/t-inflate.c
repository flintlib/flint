/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_inflate, state)
{
    slong i, j;

    /* regression test: the input used to be ignored */
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, c;
        fmpz shift[2] = {0, 0}, stride[2] = {2, 3};
        const char * vars[] = {"x", "y"};

        fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(c, ctx);

        fmpq_mpoly_set_str_pretty(b, "x + y/2", vars, ctx);
        fmpq_mpoly_set_str_pretty(c, "x^2 + y^3/2", vars, ctx);
        fmpq_mpoly_inflate(a, b, shift, stride, ctx);
        if (!fmpq_mpoly_equal(a, c, ctx))
            TEST_FUNCTION_FAIL("regression example\n");

        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(c, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* compare with substituting x_j -> x_j^stride[j] and multiplying
       by x^shift */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, m, t;
        fmpq_mpoly_struct ** C;
        fmpz * shift, * stride;
        slong nvars;

        fmpq_mpoly_ctx_init_rand(ctx, state, 5);
        nvars = ctx->zctx->minfo->nvars;
        if (nvars < 1)
        {
            fmpq_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(m, ctx);
        fmpq_mpoly_init(t, ctx);
        shift = _fmpz_vec_init(nvars);
        stride = _fmpz_vec_init(nvars);
        C = (fmpq_mpoly_struct **) flint_malloc(nvars*sizeof(fmpq_mpoly_struct *));

        fmpq_mpoly_randtest_bound(f, state, n_randint(state, 20), 50, 10, ctx);
        fmpq_mpoly_randtest_bound(g, state, n_randint(state, 20), 50, 10, ctx);

        fmpq_mpoly_one(m, ctx);
        for (j = 0; j < nvars; j++)
        {
            fmpz_set_ui(shift + j, n_randint(state, 4));
            fmpz_set_ui(stride + j, n_randint(state, 4));

            C[j] = (fmpq_mpoly_struct *) flint_malloc(sizeof(fmpq_mpoly_struct));
            fmpq_mpoly_init(C[j], ctx);
            fmpq_mpoly_gen(t, j, ctx);
            fmpq_mpoly_pow_ui(C[j], t, fmpz_get_ui(stride + j), ctx);
            fmpq_mpoly_pow_ui(t, t, fmpz_get_ui(shift + j), ctx);
            fmpq_mpoly_mul(m, m, t, ctx);
        }

        if (!fmpq_mpoly_compose_fmpq_mpoly(h, f, C, ctx, ctx))
            TEST_FUNCTION_FAIL("compose\n");
        fmpq_mpoly_mul(h, h, m, ctx);

        fmpq_mpoly_inflate(g, f, shift, stride, ctx);
        fmpq_mpoly_assert_canonical(g, ctx);
        if (!fmpq_mpoly_equal(g, h, ctx))
            TEST_FUNCTION_FAIL("i = %wd\n", i);

        fmpq_mpoly_inflate(f, f, shift, stride, ctx);
        fmpq_mpoly_assert_canonical(f, ctx);
        if (!fmpq_mpoly_equal(f, h, ctx))
            TEST_FUNCTION_FAIL("aliasing, i = %wd\n", i);

        for (j = 0; j < nvars; j++)
        {
            fmpq_mpoly_clear(C[j], ctx);
            flint_free(C[j]);
        }
        flint_free(C);
        _fmpz_vec_clear(shift, nvars);
        _fmpz_vec_clear(stride, nvars);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
