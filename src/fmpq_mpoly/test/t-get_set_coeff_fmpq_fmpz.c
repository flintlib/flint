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

TEST_FUNCTION_START(fmpq_mpoly_get_set_coeff_fmpq_fmpz, state)
{
    int i, j, k, result;

    /* Check _fmpq_fmpz */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        flint_bitcnt_t coeff_bits, exp_bits;
        slong len;
        fmpq_t c, d;

        fmpq_init(c);
        fmpq_init(d);

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);

        len = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);
        exp_bits = n_randint(state, 200) + 2;

        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz * exp = (fmpz *) flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz));

            fmpq_randtest(c, state, 200);
            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
            {
                fmpz_init(exp + k);
                fmpz_randtest_unsigned(exp + k, state, 200);
            }

            _fmpq_mpoly_set_coeff_fmpq_fmpz(f, c, exp, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);
            _fmpq_mpoly_get_coeff_fmpq_fmpz(d, f, exp, ctx);
            result = fmpq_equal(c, d);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check _fmpq_fmpz\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
                fmpz_clear(exp + k);

            flint_free(exp);
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_clear(c);
        fmpq_clear(d);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
