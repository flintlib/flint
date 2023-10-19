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

TEST_FUNCTION_START(fmpq_mpoly_get_set_coeff_fmpq_monomial, state)
{
    slong i, j, k;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, m;
        fmpz ** exp;
        fmpq_t cm, ce, q;
        slong nvars, len;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(m, ctx);
        fmpq_init(cm);
        fmpq_init(ce);
        fmpq_init(q);
        fmpq_one(q); /* anything nonzero is ok */

        nvars = fmpq_mpoly_ctx_nvars(ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, FLINT_BITS + 10) + 1;
        exp_bits1 = n_randint(state, 200);
        exp_bits2 = n_randint(state, 200);
        coeff_bits = n_randint(state, 200);
        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        /* check a random monomial - this also randomizes m->bits */
        exp = (fmpz **) flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz *));
        for (k = 0; k < nvars; k++)
        {
            exp[k] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(exp[k]);
            fmpz_randtest_unsigned(exp[k], state, exp_bits1);
        }
        fmpq_mpoly_zero(m, ctx);
        fmpq_mpoly_set_coeff_fmpq_fmpz(m, q, exp, ctx);
        fmpq_mpoly_get_coeff_fmpq_monomial(cm, f, m, ctx);
        fmpq_mpoly_get_coeff_fmpq_fmpz(ce, f, exp, ctx);
        if (!fmpq_equal(cm, ce))
        {
            flint_printf("FAIL\ncheck a random monomial\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        /* check all monomials in f */
        for (j = 0; j < fmpq_mpoly_length(f, ctx); j++)
        {
            fmpq_mpoly_get_term_exp_fmpz(exp, f, j, ctx);

            fmpq_mpoly_zero(m, ctx);
            fmpq_mpoly_set_coeff_fmpq_fmpz(m, q, exp, ctx);
            fmpq_mpoly_get_coeff_fmpq_monomial(cm, f, m, ctx);
            fmpq_mpoly_get_coeff_fmpq_fmpz(ce, f, exp, ctx);
            if (!fmpq_equal(cm, ce))
            {
                flint_printf("FAIL\ncheck all monomials in f\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        /* set random coeff and check */
        for (j = 0; j < 10; j++)
        {
            for (k = 0; k < nvars; k++)
            {
                fmpz_randtest_unsigned(exp[k], state, exp_bits2);
            }

            fmpq_randtest(cm, state, coeff_bits + 1);

            fmpq_mpoly_zero(m, ctx);
            fmpq_mpoly_set_coeff_fmpq_fmpz(m, q, exp, ctx);
            fmpq_mpoly_set_coeff_fmpq_monomial(f, cm, m, ctx);
            fmpq_mpoly_get_coeff_fmpq_monomial(ce, f, m, ctx);
            if (!fmpq_equal(cm, ce))
            {
                flint_printf("FAIL\nset random coeff and check\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (k = 0; k < nvars; k++)
        {
            fmpz_clear(exp[k]);
            flint_free(exp[k]);
        }
        flint_free(exp);

        fmpq_clear(q);
        fmpq_clear(cm);
        fmpq_clear(ce);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
