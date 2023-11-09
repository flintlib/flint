/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_degree, state)
{
    int i, j;

    /* Check degree does not go up under addition */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 4);

        for (j = 0; j < ctx->zctx->minfo->nvars; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_add(h, f, g, ctx);

            fmpq_mpoly_degree_fmpz(hdeg, h, j, ctx);
            fmpq_mpoly_degree_fmpz(fdeg, f, j, ctx);
            fmpq_mpoly_degree_fmpz(gdeg, g, j, ctx);

            if ((fmpz_cmp(hdeg, fdeg) > 0) && (fmpz_cmp(hdeg, gdeg) > 0))
            {
                printf("FAIL\n");
                flint_printf("Check degree does not go up under addition\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check degree adds under multiplication */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        int ok;
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 4);

        for (j = 0; j < ctx->zctx->minfo->nvars; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_mul(h, f, g, ctx);

            fmpq_mpoly_degree_fmpz(hdeg, h, j, ctx);
            fmpq_mpoly_degree_fmpz(fdeg, f, j, ctx);
            fmpq_mpoly_degree_fmpz(gdeg, g, j, ctx);

            if (fmpq_mpoly_is_zero(g, ctx) || fmpq_mpoly_is_zero(f, ctx)) {
                ok = (fmpz_cmp_si(hdeg, -WORD(1)) == 0);
            } else {
                fmpz_sub(hdeg, hdeg, fdeg);
                ok = (fmpz_cmp(hdeg, gdeg) == 0);
            }

            if (!ok)
            {
                printf("FAIL\n");
                flint_printf("Check degree adds under multiplication\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
