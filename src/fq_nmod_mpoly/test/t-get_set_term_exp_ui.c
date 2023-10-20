/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_get_set_term_exp_ui, state)
{
    slong i, j, k;
    int result;

    /* check set and get match */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f;
        slong nvars, len, index;
        flint_bitcnt_t exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);

        nvars = fq_nmod_mpoly_ctx_nvars(ctx);

        len = n_randint(state, 50) + 1;
        exp_bits = n_randint(state, 100) + 2;

        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
        if (fq_nmod_mpoly_is_zero(f, ctx))
            fq_nmod_mpoly_one(f, ctx);

        for (j = 0; j < 10; j++)
        {
            ulong * exp1 = (ulong *) flint_malloc(nvars*sizeof(ulong));
            ulong * exp2 = (ulong *) flint_malloc(nvars*sizeof(ulong));

            for (k = 0; k < nvars; k++)
            {
                slong bits = n_randint(state, FLINT_BITS) + 1;
                exp1[k] = n_randbits(state, bits);
            }

            index = n_randint(state, f->length);

            fq_nmod_mpoly_set_term_exp_ui(f, index, exp1, ctx);

            if (!mpoly_monomials_valid_test(f->exps, f->length, f->bits, ctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

            if (mpoly_monomials_overflow_test(f->exps, f->length, f->bits, ctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

            fq_nmod_mpoly_get_term_exp_ui(exp2, f, index, ctx);

            result = 1;
            for (k = 0; k < nvars; k++)
            {
                result = result
                 && exp1[k] == exp2[k]
                 && exp1[k] == fq_nmod_mpoly_get_term_var_exp_ui(f, index, k, ctx);
            }

            if (!result)
            {
                flint_printf("FAIL\ncheck set and get match\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            flint_free(exp1);
            flint_free(exp2);
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
