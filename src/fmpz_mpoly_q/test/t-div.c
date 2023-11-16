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

TEST_FUNCTION_START(fmpz_mpoly_q_div, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_t A, B, C, D;

        fmpz_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX);

        fmpz_mpoly_q_init(A, ctx);
        fmpz_mpoly_q_init(B, ctx);
        fmpz_mpoly_q_init(C, ctx);
        fmpz_mpoly_q_init(D, ctx);

        fmpz_mpoly_q_randtest(A, state, 10, 2 + n_randint(state, 100), 5, ctx);
        do {
            fmpz_mpoly_q_randtest(B, state, 10, 2 + n_randint(state, 100), 5, ctx);
        } while (fmpz_mpoly_q_is_zero(B, ctx));

        fmpz_mpoly_q_div(C, A, B, ctx);

        fmpz_mpoly_mul(fmpz_mpoly_q_numref(D), fmpz_mpoly_q_numref(A), fmpz_mpoly_q_denref(B), ctx);
        fmpz_mpoly_mul(fmpz_mpoly_q_denref(D), fmpz_mpoly_q_denref(A), fmpz_mpoly_q_numref(B), ctx);
        fmpz_mpoly_q_canonicalise(D, ctx);

        if (!fmpz_mpoly_q_equal(C, D, ctx))
        {
            flint_printf("FAIL\n");
            flint_printf("A = "); fmpz_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("D = "); fmpz_mpoly_q_print_pretty(D, NULL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            fmpz_mpoly_q_set(C, A, ctx);
            fmpz_mpoly_q_div(C, C, B, ctx);
        }
        else
        {
            fmpz_mpoly_q_set(C, B, ctx);
            fmpz_mpoly_q_div(C, A, C, ctx);
        }

        if (!fmpz_mpoly_q_equal(C, D, ctx))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("A = "); fmpz_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("D = "); fmpz_mpoly_q_print_pretty(D, NULL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mpoly_q_clear(A, ctx);
        fmpz_mpoly_q_clear(B, ctx);
        fmpz_mpoly_q_clear(C, ctx);
        fmpz_mpoly_q_clear(D, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
