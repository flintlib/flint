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
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpq_mpoly_get_set_term_exp_si, state)
{
    slong i, j, k;
    int result;

    /* check get and set match */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        slong len, nvars, index;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);

        nvars = fmpq_mpoly_ctx_nvars(ctx);

        len = n_randint(state, 50) + 1;
        exp_bits = n_randint(state, 200) + 2;
        coeff_bits = n_randint(state, 100) + 1;

        do {
            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
        } while (fmpq_mpoly_length(f, ctx) == 0);

        for (j = 0; j < 10; j++)
        {
            slong * exp1 = (slong *) flint_malloc(ctx->zctx->minfo->nvars*sizeof(slong));
            slong * exp2 = (slong *) flint_malloc(ctx->zctx->minfo->nvars*sizeof(slong));

            for (k = 0; k < nvars; k++)
            {
                slong bits = n_randint(state, FLINT_BITS - 1) + 1;
                exp1[k] = n_randbits(state, bits);
                exp1[k] = FLINT_ABS(exp1[k]);
            }

            index = n_randint(state, fmpq_mpoly_length(f, ctx));

            fmpq_mpoly_set_term_exp_ui(f, index, (ulong *) exp1, ctx);

            if (!mpoly_monomials_valid_test(f->zpoly->exps, f->zpoly->length, f->zpoly->bits, ctx->zctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

            if (mpoly_monomials_overflow_test(f->zpoly->exps, f->zpoly->length, f->zpoly->bits, ctx->zctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

            fmpq_mpoly_get_term_exp_si(exp2, f, index, ctx);

            result = 1;
            for (k = 0; k < nvars; k++)
            {
                result = result
                 && exp1[k] == exp2[k]
                 && exp1[k] == fmpq_mpoly_get_term_var_exp_si(f, index, k, ctx);
            }

            if (!result)
            {
                flint_printf("FAIL\ncheck get and set match\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            flint_free(exp1);
            flint_free(exp2);
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
