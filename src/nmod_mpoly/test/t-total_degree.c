/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

/* foolproof way to check totdeg_check is correct */
void _check_total_degree(const fmpz_t totdeg_check, const nmod_mpoly_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    fmpz_t totdeg;
    fmpz_init(totdeg);
    mpoly_total_degree_fmpz_ref(totdeg, A->exps, A->length, A->bits, ctx->minfo);
    if (!fmpz_equal(totdeg_check, totdeg))
        flint_throw(FLINT_ERROR, "Total degree is wrong");
    fmpz_clear(totdeg);
}

TEST_FUNCTION_START(nmod_mpoly_total_degree, state)
{
    int i, j;

    /* Check total_degree does not go up under addition */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randbits(state, n_randint(state, FLINT_BITS));
        modulus = FLINT_MAX(UWORD(2), modulus);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_add(h, f, g, ctx);

            nmod_mpoly_total_degree_fmpz(hdeg, h, ctx);
            nmod_mpoly_total_degree_fmpz(fdeg, f, ctx);
            nmod_mpoly_total_degree_fmpz(gdeg, g, ctx);
            _check_total_degree(hdeg, h, ctx);
            _check_total_degree(fdeg, f, ctx);
            _check_total_degree(gdeg, g, ctx);

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
        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check total_degree adds under multiplication */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        int ok;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randbits(state, n_randint(state, FLINT_BITS));
        modulus = FLINT_MAX(UWORD(2), modulus);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_mul_johnson(h, f, g, ctx);

            nmod_mpoly_total_degree_fmpz(hdeg, h, ctx);
            nmod_mpoly_total_degree_fmpz(fdeg, f, ctx);
            nmod_mpoly_total_degree_fmpz(gdeg, g, ctx);
            _check_total_degree(hdeg, h, ctx);
            _check_total_degree(fdeg, f, ctx);
            _check_total_degree(gdeg, g, ctx);

            if (nmod_mpoly_is_zero(f, ctx) || nmod_mpoly_is_zero(g, ctx))
            {
                ok = fmpz_equal_si(hdeg, -WORD(1))
                     && (fmpz_equal_si(fdeg, -WORD(1))
                         || fmpz_equal_si(gdeg, -WORD(1)));
            }
            else
            {
                fmpz_add(gdeg, gdeg, fdeg);
                if (n_is_prime(nmod_mpoly_ctx_modulus(ctx)))
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
                printf("FAIL\n");
                flint_printf("Check degree adds under multiplication\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
