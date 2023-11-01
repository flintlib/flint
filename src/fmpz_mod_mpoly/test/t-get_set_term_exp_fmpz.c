/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_get_set_term_exp_fmpz, state)
{
    slong i, j, k;
    int result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f;
        slong nvars, len, index;
        flint_bitcnt_t exp_bits;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        fmpz_mod_mpoly_init(f, ctx);

        nvars = fmpz_mod_mpoly_ctx_nvars(ctx);

        do {
            len = n_randint(state, 50);
            exp_bits = n_randint(state, 100) + 1;

            fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
        } while (f->length == 0);

        for (j = 0; j < 10; j++)
        {
            fmpz ** exp1 = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
            fmpz ** exp2 = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));

            for (k = 0; k < nvars; k++)
            {
                exp1[k] = (fmpz *) flint_malloc(sizeof(fmpz));
                exp2[k] = (fmpz *) flint_malloc(sizeof(fmpz));
                fmpz_init(exp1[k]);
                fmpz_init(exp2[k]);
                fmpz_randtest_unsigned(exp1[k], state, 200);
            }

            index = n_randint(state, f->length);

            fmpz_mod_mpoly_set_term_exp_fmpz(f, index, exp1, ctx);

            if (!mpoly_monomials_valid_test(f->exps, f->length, f->bits, ctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

            if (mpoly_monomials_overflow_test(f->exps, f->length, f->bits, ctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

            fmpz_mod_mpoly_get_term_exp_fmpz(exp2, f, index, ctx);

            result = 1;
            for (k = 0; k < nvars; k++)
            {
                result = result && fmpz_equal(exp1[k], exp2[k]);
                fmpz_clear(exp1[k]);
                fmpz_clear(exp2[k]);
                flint_free(exp1[k]);
                flint_free(exp2[k]);
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

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
