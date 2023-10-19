/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_mpoly_get_set_coeff_fmpz_monomial, state)
{
    slong i, j, k;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, m;
        fmpz ** exp;
        fmpz_t cm, ce, q;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(m, ctx);
        fmpz_init(cm);
        fmpz_init(ce);
        fmpz_init(q);
        fmpz_one(q); /* anything nonzero is ok */

        len = n_randint(state, 100);
        exp_bits = n_randint(state, FLINT_BITS + 10) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200);
        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        /* check a random monomial - this also randomizes m->bits */
        exp = (fmpz **) flint_malloc(ctx->minfo->nvars*sizeof(fmpz *));
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            exp[k] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(exp[k]);
            fmpz_randtest_unsigned(exp[k], state, exp_bits1);
        }
        fmpz_mpoly_zero(m, ctx);
        fmpz_mpoly_set_coeff_fmpz_fmpz(m, q, exp, ctx);
        fmpz_mpoly_get_coeff_fmpz_monomial(cm, f, m, ctx);
        fmpz_mpoly_get_coeff_fmpz_fmpz(ce, f, exp, ctx);
        if (!fmpz_equal(cm, ce))
        {
            flint_printf("FAIL\ncheck a random monomial\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        /* check all coeffs in f */
        for (j = 0; j < fmpz_mpoly_length(f, ctx); j++)
        {
            fmpz_mpoly_get_term_exp_fmpz(exp, f, j, ctx);

            fmpz_mpoly_zero(m, ctx);
            fmpz_mpoly_set_coeff_fmpz_fmpz(m, q, exp, ctx);
            fmpz_mpoly_get_coeff_fmpz_monomial(cm, f, m, ctx);
            fmpz_mpoly_get_coeff_fmpz_fmpz(ce, f, exp, ctx);
            if (!fmpz_equal(cm, ce))
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

            fmpz_randtest(cm, state, coeff_bits);

            fmpz_mpoly_zero(m, ctx);
            fmpz_mpoly_set_coeff_fmpz_fmpz(m, q, exp, ctx);
            fmpz_mpoly_set_coeff_fmpz_monomial(f, cm, m, ctx);
            fmpz_mpoly_get_coeff_fmpz_monomial(ce, f, m, ctx);
            if (!fmpz_equal(cm, ce))
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

        fmpz_clear(q);
        fmpz_clear(cm);
        fmpz_clear(ce);
        fmpz_mpoly_clear(m, ctx);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
