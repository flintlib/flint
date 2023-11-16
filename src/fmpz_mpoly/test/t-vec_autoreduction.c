/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_vec_autoreduction, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_vec_t F, G;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);

        fmpz_mpoly_vec_init(F, 0, ctx);
        fmpz_mpoly_vec_init(G, 0, ctx);

        fmpz_mpoly_vec_randtest_not_zero(F, state, 1 + n_randint(state, 5), 1 + n_randint(state, 4), 1 + n_randint(state, 3), 1 + n_randint(state, 2), ctx);

        fmpz_mpoly_vec_autoreduction(G, F, ctx);

        if (!fmpz_mpoly_vec_is_autoreduced(G, ctx))
        {
            flint_printf("FAIL\n\n");
            mpoly_ordering_print(ctx->minfo->ord); printf("\n");
            flint_printf("F = "); fmpz_mpoly_vec_print(F, ctx); flint_printf("\n");
            flint_printf("G = "); fmpz_mpoly_vec_print(G, ctx); flint_printf("\n");
            flint_abort();
        }

        fmpz_mpoly_vec_clear(F, ctx);
        fmpz_mpoly_vec_clear(G, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
