/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_pow_ui, state)
{
    int i, j, k, result;

    /* Check against rmul */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len1, len2;
        fmpz_t power;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpz_init(power);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 10; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            fmpq_mpoly_pow_ui(h, f, j, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            fmpq_mpoly_one(g, ctx);
            for (k = 0; k < j; k++)
                fmpq_mpoly_mul(g, g, f, ctx);

            result = fmpq_mpoly_equal(h, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check against rmul\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(power);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check monomials against pow_fmpz */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len2;
        fmpz_t power;
        flint_bitcnt_t coeff_bits, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        fmpz_init(power);

        len2 = n_randint(state, 10);
        exp_bits2 = n_randint(state, 200) + 2;
        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 10; j++)
        {
            /* make sure power is random ui */
            fmpz_set_ui(power, n_randlimb(state));

            /* set f to a random monomial */
            fmpq_mpoly_one(f, ctx);
            if (n_randint(state, 2))
            {
                fmpq_mpoly_neg(f, f, ctx);
            }
            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
            {
                fmpq_mpoly_gen(h, n_randint(state, ctx->zctx->minfo->nvars), ctx);
                fmpq_mpoly_mul(f, f, h, ctx);
            }

            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            fmpq_mpoly_pow_ui(g, f, fmpz_get_ui(power), ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            fmpq_mpoly_pow_fmpz(h, f, power, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            result = fmpq_mpoly_equal(h, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check monomials against pow_fmpz\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(power);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
