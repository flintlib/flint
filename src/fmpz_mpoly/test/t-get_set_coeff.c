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
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_get_set_coeff, state)
{
    int i, j, k, result;

    /* Check _fmpz_fmpz */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits2;
        slong len;
        fmpz_t c, d;

        fmpz_init(c);
        fmpz_init(d);

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);

        len = n_randint(state, 100);
        coeff_bits = n_randint(state, 200);
        exp_bits = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz ** exp = (fmpz **) flint_malloc(ctx->minfo->nvars*sizeof(fmpz*));

            fmpz_randtest_unsigned(c, state, 200);
            for (k = 0; k < fmpz_mpoly_ctx_nvars(ctx); k++)
            {
                exp[k] = (fmpz *) flint_malloc(sizeof(fmpz));
                fmpz_init(exp[k]);
                fmpz_randtest_unsigned(exp[k], state, exp_bits2);
            }

            fmpz_mpoly_set_coeff_fmpz_fmpz(f, c, exp, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            fmpz_mpoly_get_coeff_fmpz_fmpz(d, f, exp, ctx);
            result = fmpz_equal(c, d);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check _fmpz_fmpz\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            for (k = 0; k < fmpz_mpoly_ctx_nvars(ctx); k++)
            {
                fmpz_clear(exp[k]);
                flint_free(exp[k]);
            }

            flint_free(exp);
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check _fmpz_ui */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        flint_bitcnt_t coeff_bits, exp_bits;
        slong len;
        fmpz_t c, d;

        fmpz_init(c);
        fmpz_init(d);

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);

        len = n_randint(state, 100);
        coeff_bits = n_randint(state, 200);
        exp_bits = n_randint(state, 200) + 2;

        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            ulong * exp = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

            fmpz_randtest_unsigned(c, state, 200);
            for (k = 0; k < fmpz_mpoly_ctx_nvars(ctx); k++)
                exp[k] = n_randtest(state);

            fmpz_mpoly_set_coeff_fmpz_ui(f, c, exp, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            fmpz_mpoly_get_coeff_fmpz_ui(d, f, exp, ctx);
            result = fmpz_equal(c, d);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check _fmpz_ui\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            flint_free(exp);
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
