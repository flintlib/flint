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
#include "fmpz_poly.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_get_set_is_fmpz_poly, state)
{
    slong i;

    /* regression test: inputs that are not univariate */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a;
        fmpz_poly_t p;
        const char * vars[] = {"x", "y"};

        fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_poly_init(p);

        fmpz_mpoly_set_str_pretty(a, "x*y + 2*x + 3*y + 5", vars, ctx);
        if (fmpz_mpoly_get_fmpz_poly(p, a, 0, ctx) ||
            fmpz_mpoly_get_fmpz_poly(p, a, 1, ctx))
        {
            TEST_FUNCTION_FAIL("not univariate\n");
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_poly_clear(p);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b;
        fmpz_poly_t p, q;
        slong var;
        int r;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        if (ctx->minfo->nvars < 1)
        {
            fmpz_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_poly_init(p);
        fmpz_poly_init(q);

        var = n_randint(state, ctx->minfo->nvars);

        if (n_randint(state, 2))
        {
            /* a univariate polynomial, possibly with multiword exponents */
            fmpz_poly_randtest(p, state, n_randint(state, 20), 50);
            fmpz_mpoly_set_fmpz_poly(a, p, var, ctx);
            if (n_randint(state, 2))
            {
                flint_bitcnt_t bits = mpoly_fix_bits(FLINT_BITS + 1 +
                                     n_randint(state, FLINT_BITS), ctx->minfo);
                if (!fmpz_mpoly_repack_bits(a, a, bits, ctx))
                    TEST_FUNCTION_FAIL("repack_bits\n");
            }

            if (!fmpz_mpoly_is_fmpz_poly(a, var, ctx) ||
                !fmpz_mpoly_get_fmpz_poly(q, a, var, ctx) ||
                !fmpz_poly_equal(p, q))
            {
                TEST_FUNCTION_FAIL("univariate input, i = %wd\n", i);
            }
        }
        else
        {
            /* small exponents, possibly stored with multiword exponents */
            fmpz_mpoly_randtest_bits(a, state, n_randint(state, 5), 50,
                                                 1 + n_randint(state, 3), ctx);
            if (n_randint(state, 2))
            {
                flint_bitcnt_t bits = mpoly_fix_bits(FLINT_BITS + 1 +
                               n_randint(state, FLINT_BITS), ctx->minfo);
                if (!fmpz_mpoly_repack_bits(a, a, bits, ctx))
                    TEST_FUNCTION_FAIL("repack_bits\n");
            }

            r = fmpz_mpoly_get_fmpz_poly(q, a, var, ctx);
            if (r != fmpz_mpoly_is_fmpz_poly(a, var, ctx))
            {
                TEST_FUNCTION_FAIL("return value, i = %wd\n", i);
            }

            if (r)
            {
                fmpz_mpoly_set_fmpz_poly(b, q, var, ctx);
                if (!fmpz_mpoly_equal(a, b, ctx))
                    TEST_FUNCTION_FAIL("get/set, i = %wd\n", i);
            }
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
