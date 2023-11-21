/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly_q.h"

TEST_FUNCTION_START(fmpz_mpoly_q_randtest, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_t A;

        fmpz_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX);

        fmpz_mpoly_q_init(A, ctx);
        fmpz_mpoly_q_randtest(A, state, 10, 10, 5, ctx);

        /*
            fmpz_mpoly_q_print_pretty(A, NULL, ctx);
            flint_printf("\n");
        */

        if (!fmpz_mpoly_q_is_canonical(A, ctx))
        {
            flint_printf("FAIL: not canonical\n");
            fmpz_mpoly_q_print_pretty(A, NULL, ctx);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mpoly_q_clear(A, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
