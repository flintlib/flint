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

TEST_FUNCTION_START(fmpz_mod_mpoly_get_set_coeff_fmpz_fmpz, state)
{
    slong i, j, k;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f;
        flint_bitcnt_t exp_bits, exp_bits2;
        slong len;
        fmpz_t c, d;

        fmpz_init(c);
        fmpz_init(d);

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        fmpz_mod_mpoly_init(f, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz ** exp = (fmpz **) flint_malloc(ctx->minfo->nvars*sizeof(fmpz*));

            fmpz_randtest_unsigned(c, state, 200);
            for (k = 0; k < fmpz_mod_mpoly_ctx_nvars(ctx); k++)
            {
                exp[k] = (fmpz *) flint_malloc(sizeof(fmpz));
                fmpz_init(exp[k]);
                fmpz_randtest_unsigned(exp[k], state, exp_bits2);
            }

            fmpz_mod_mpoly_set_coeff_fmpz_fmpz(f, c, exp, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            fmpz_mod_mpoly_get_coeff_fmpz_fmpz(d, f, exp, ctx);

            if (!fmpz_mod_equal_fmpz(d, c, ctx->ffinfo))
            {
                printf("FAIL\n");
                flint_printf("Check _fmpz_fmpz\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            for (k = 0; k < fmpz_mod_mpoly_ctx_nvars(ctx); k++)
            {
                fmpz_clear(exp[k]);
                flint_free(exp[k]);
            }

            flint_free(exp);
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
