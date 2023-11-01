/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_mul_monomial, state)
{
    slong i, j;
    slong tmul = 200;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, d, bb, cc;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(d, ctx);
        fmpz_mpoly_init(bb, ctx);
        fmpz_mpoly_init(cc, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz_mpoly_randtest_bits(a, state, 1 + n_randint(state, 50),
                                               1 + n_randint(state, 200),
                                               1 + n_randint(state, 192), ctx);

            fmpz_mpoly_randtest_bits(b, state, 1 + n_randint(state, 30),
                                               1 + n_randint(state, 200),
                                               1 + n_randint(state, 192), ctx);

            fmpz_mpoly_randtest_bits(c, state, 1,
                                               1 + n_randint(state, 200),
                                               1 + n_randint(state, 192), ctx);

            if (fmpz_mpoly_length(c, ctx) != 1 || n_randint(state, 50) == 0)
                fmpz_mpoly_set_ui(c, n_randtest_not_zero(state), ctx);

            fmpz_mpoly_mul_monomial(a, b, c, ctx);
            fmpz_mpoly_assert_canonical(a, ctx);
            fmpz_mpoly_mul_johnson(d, b, c, ctx);
            if (!fmpz_mpoly_equal(a, d, ctx))
            {
                flint_printf("FAIL: check mul_monomial against mul_johnson\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mpoly_set(bb, b, ctx);
            fmpz_mpoly_set(cc, c, ctx);
            fmpz_mpoly_mul_monomial(bb, bb, cc, ctx);
            fmpz_mpoly_assert_canonical(bb, ctx);
            if (!fmpz_mpoly_equal(bb, d, ctx) ||
                 !fmpz_mpoly_equal(cc, c, ctx))
            {
                flint_printf("FAIL: check aliasing first input\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mpoly_set(bb, b, ctx);
            fmpz_mpoly_set(cc, c, ctx);
            fmpz_mpoly_mul_monomial(cc, bb, cc, ctx);
            fmpz_mpoly_assert_canonical(cc, ctx);
            if (!fmpz_mpoly_equal(cc, d, ctx) ||
                 !fmpz_mpoly_equal(bb, b, ctx))
            {
                flint_printf("FAIL: check aliasing second input\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mpoly_set(cc, c, ctx);
            fmpz_mpoly_mul_monomial(cc, cc, cc, ctx);
            fmpz_mpoly_assert_canonical(cc, ctx);
            fmpz_mpoly_mul_johnson(d, c, c, ctx);
            if (!fmpz_mpoly_equal(cc, d, ctx))
            {
                flint_printf("FAIL: check aliasing both inputs\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(d, ctx);
        fmpz_mpoly_clear(bb, ctx);
        fmpz_mpoly_clear(cc, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
