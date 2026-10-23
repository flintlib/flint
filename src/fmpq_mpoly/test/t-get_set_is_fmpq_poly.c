/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpoly.h"
#include "fmpq_poly.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_get_set_is_fmpq_poly, state)
{
    slong i;

    /* regression test: inputs that are not univariate */
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a;
        fmpq_poly_t p;
        const char * vars[] = {"x", "y"};

        fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
        fmpq_mpoly_init(a, ctx);
        fmpq_poly_init(p);

        fmpq_mpoly_set_str_pretty(a, "x*y + 2*x + 3*y + 5", vars, ctx);
        if (fmpq_mpoly_get_fmpq_poly(p, a, 0, ctx) ||
            fmpq_mpoly_get_fmpq_poly(p, a, 1, ctx))
        {
            TEST_FUNCTION_FAIL("not univariate\n");
        }

        fmpq_mpoly_clear(a, ctx);
        fmpq_poly_clear(p);
        fmpq_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b;
        fmpq_poly_t p, q;
        slong var;
        int r;

        fmpq_mpoly_ctx_init_rand(ctx, state, 4);
        if (ctx->zctx->minfo->nvars < 1)
        {
            fmpq_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_poly_init(p);
        fmpq_poly_init(q);

        var = n_randint(state, ctx->zctx->minfo->nvars);

        if (n_randint(state, 2))
        {
            /* a univariate polynomial, possibly with multiword exponents */
            fmpq_poly_randtest(p, state, n_randint(state, 20), 50);
            fmpq_mpoly_set_fmpq_poly(a, p, var, ctx);
            if (n_randint(state, 2))
            {
                flint_bitcnt_t bits = mpoly_fix_bits(FLINT_BITS + 1 +
                               n_randint(state, FLINT_BITS), ctx->zctx->minfo);
                if (!fmpq_mpoly_repack_bits(a, a, bits, ctx))
                    TEST_FUNCTION_FAIL("repack_bits\n");
            }

            if (!fmpq_mpoly_is_fmpq_poly(a, var, ctx) ||
                !fmpq_mpoly_get_fmpq_poly(q, a, var, ctx) ||
                !fmpq_poly_equal(p, q))
            {
                TEST_FUNCTION_FAIL("univariate input, i = %wd\n", i);
            }
        }
        else
        {
            /* small exponents, possibly stored with multiword exponents */
            fmpq_mpoly_randtest_bits(a, state, n_randint(state, 5), 50,
                                                 1 + n_randint(state, 3), ctx);
            if (n_randint(state, 2))
            {
                flint_bitcnt_t bits = mpoly_fix_bits(FLINT_BITS + 1 +
                               n_randint(state, FLINT_BITS), ctx->zctx->minfo);
                if (!fmpq_mpoly_repack_bits(a, a, bits, ctx))
                    TEST_FUNCTION_FAIL("repack_bits\n");
            }

            r = fmpq_mpoly_get_fmpq_poly(q, a, var, ctx);
            if (r != fmpq_mpoly_is_fmpq_poly(a, var, ctx))
            {
                TEST_FUNCTION_FAIL("return value, i = %wd\n", i);
            }

            if (r)
            {
                fmpq_mpoly_set_fmpq_poly(b, q, var, ctx);
                if (!fmpq_mpoly_equal(a, b, ctx))
                    TEST_FUNCTION_FAIL("get/set, i = %wd\n", i);
            }
        }

        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_poly_clear(p);
        fmpq_poly_clear(q);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
