/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_get_set_coeff, state)
{
    int i, j, k, result;

    /* Check _fmpz_fmpz */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f;
        flint_bitcnt_t exp_bits, exp_bits2;
        slong len;
        fq_nmod_t c, d;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_init(d, ctx->fqctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz ** exp = (fmpz **) flint_malloc(ctx->minfo->nvars*sizeof(fmpz*));

            fq_nmod_randtest(c, state, ctx->fqctx);
            for (k = 0; k < fq_nmod_mpoly_ctx_nvars(ctx); k++)
            {
                exp[k] = (fmpz *) flint_malloc(sizeof(fmpz));
                fmpz_init(exp[k]);
                fmpz_randtest_unsigned(exp[k], state, exp_bits2);
            }

            fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(f, c, exp, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);
            fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(d, f, exp, ctx);
            result = fq_nmod_equal(c, d, ctx->fqctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check _fmpz_fmpz\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            for (k = 0; k < fq_nmod_mpoly_ctx_nvars(ctx); k++)
            {
                fmpz_clear(exp[k]);
                flint_free(exp[k]);
            }

            flint_free(exp);
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_clear(d, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check _fmpz_ui */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f;
        flint_bitcnt_t exp_bits;
        slong len;
        fq_nmod_t c, d;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_init(d, ctx->fqctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;

        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            ulong * exp = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

            fq_nmod_randtest(c, state, ctx->fqctx);
            for (k = 0; k < fq_nmod_mpoly_ctx_nvars(ctx); k++)
                exp[k] = n_randtest(state);

            fq_nmod_mpoly_set_coeff_fq_nmod_ui(f, c, exp, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);
            fq_nmod_mpoly_get_coeff_fq_nmod_ui(d, f, exp, ctx);
            result = fq_nmod_equal(c, d, ctx->fqctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check _fmpz_ui\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            flint_free(exp);
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_clear(d, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
