/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_degree, state)
{
    slong i, j;

    /* Check degree does not go up under addition */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 100);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fmpz_mod_mpoly_add(h, f, g, ctx);

            fmpz_mod_mpoly_degree_fmpz(hdeg, h, j, ctx);
            fmpz_mod_mpoly_degree_fmpz(fdeg, f, j, ctx);
            fmpz_mod_mpoly_degree_fmpz(gdeg, g, j, ctx);

            if ((fmpz_cmp(hdeg, fdeg) > 0) && (fmpz_cmp(hdeg, gdeg) > 0))
            {
                flint_printf("FAIL: Check degree does not go up under addition\n");
                flint_printf("i: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check degree adds under multiplication */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        int ok;
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 100);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fmpz_mod_mpoly_mul_johnson(h, f, g, ctx);

            fmpz_mod_mpoly_degree_fmpz(hdeg, h, j, ctx);
            fmpz_mod_mpoly_degree_fmpz(fdeg, f, j, ctx);
            fmpz_mod_mpoly_degree_fmpz(gdeg, g, j, ctx);

            if (fmpz_mod_mpoly_is_zero(f, ctx) || fmpz_mod_mpoly_is_zero(g, ctx))
            {
                ok = fmpz_equal_si(hdeg, -WORD(1))
                     && (fmpz_equal_si(fdeg, -WORD(1))
                         || fmpz_equal_si(gdeg, -WORD(1)));
            }
            else
            {
                fmpz_add(gdeg, gdeg, fdeg);
                if (fmpz_is_prime(fmpz_mod_mpoly_ctx_modulus(ctx)))
                {
                    ok = fmpz_equal(hdeg, gdeg);
                }
                else
                {
                    ok = fmpz_cmp(hdeg, gdeg) <= 0;
                }
            }

            if (!ok)
            {
                flint_printf("FAIL: Check degree adds under multiplication\n");
                flint_printf("i: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
