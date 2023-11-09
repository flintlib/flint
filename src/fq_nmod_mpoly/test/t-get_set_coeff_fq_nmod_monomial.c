/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_get_set_coeff_fq_nmod_monomial, state)
{
    slong i, j, k;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, m;
        fmpz ** exp;
        fq_nmod_t cm, ce, q;
        slong len;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(m, ctx);
        fq_nmod_init(cm, ctx->fqctx);
        fq_nmod_init(ce, ctx->fqctx);
        fq_nmod_init(q, ctx->fqctx);
        fq_nmod_one(q, ctx->fqctx); /* anything nonzero is ok */

        len = n_randint(state, 100);
        exp_bits = n_randint(state, FLINT_BITS + 10) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;
        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        /* check a random monomial - this also randomizes m->bits */
        exp = (fmpz **) flint_malloc(ctx->minfo->nvars*sizeof(fmpz *));
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            exp[k] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(exp[k]);
            fmpz_randtest_unsigned(exp[k], state, exp_bits1);
        }
        fq_nmod_mpoly_zero(m, ctx);
        fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(m, q, exp, ctx);
        fq_nmod_mpoly_get_coeff_fq_nmod_monomial(cm, f, m, ctx);
        fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(ce, f, exp, ctx);
        if (!fq_nmod_equal(cm, ce, ctx->fqctx))
        {
            flint_printf("FAIL\ncheck a random monomial\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        /* check all coeffs in f */
        for (j = 0; j < fq_nmod_mpoly_length(f, ctx); j++)
        {
            fq_nmod_mpoly_get_term_exp_fmpz(exp, f, j, ctx);

            fq_nmod_mpoly_zero(m, ctx);
            fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(m, q, exp, ctx);
            fq_nmod_mpoly_get_coeff_fq_nmod_monomial(cm, f, m, ctx);
            fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(ce, f, exp, ctx);
            if (!fq_nmod_equal(cm, ce, ctx->fqctx))
            {
                flint_printf("FAIL\ncheck all coeffs in f\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        /* set random coeff and check */
        for (j = 0; j < 10; j++)
        {
            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(exp[k], state, exp_bits2);
            }

            fq_nmod_randtest(cm, state, ctx->fqctx);

            fq_nmod_mpoly_zero(m, ctx);
            fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(m, q, exp, ctx);
            fq_nmod_mpoly_set_coeff_fq_nmod_monomial(f, cm, m, ctx);
            fq_nmod_mpoly_get_coeff_fq_nmod_monomial(ce, f, m, ctx);
            if (!fq_nmod_equal(cm, ce, ctx->fqctx))
            {
                flint_printf("FAIL\nset random coeff and check\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_clear(exp[k]);
            flint_free(exp[k]);
        }
        flint_free(exp);

        fq_nmod_clear(q, ctx->fqctx);
        fq_nmod_clear(cm, ctx->fqctx);
        fq_nmod_clear(ce, ctx->fqctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
