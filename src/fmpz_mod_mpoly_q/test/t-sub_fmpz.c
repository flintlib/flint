/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly_q.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_q_sub_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_q_t A, B, C;
        fmpz_t c;

        fmpz_t m;

        fmpz_init(m);
        fmpz_randtest_unsigned(m, state, n_randint(state, 2) ? 4 : n_randint(state, 100));
        fmpz_nextprime(m, m, 0);
        fmpz_mod_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX, m);

        fmpz_mod_mpoly_q_init(A, ctx);
        fmpz_mod_mpoly_q_init(B, ctx);
        fmpz_mod_mpoly_q_init(C, ctx);
        fmpz_init(c);

        fmpz_mod_mpoly_q_randtest(A, state, 10, 5, ctx);
        fmpz_randtest(c, state, 10);

        fmpz_mod_mpoly_q_sub_fmpz(B, A, c, ctx);

        fmpz_mod_mpoly_q_set_fmpz(C, c, ctx);
        fmpz_mod_mpoly_q_sub(C, A, C, ctx);

        if (!fmpz_mod_mpoly_q_equal(B, C, ctx))
        {
            flint_printf("FAIL\n");
            flint_printf("A = "); fmpz_mod_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mod_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mod_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n\n");
            flint_printf("mod = "); fmpz_print(ctx->ffinfo->n); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mod_mpoly_q_clear(A, ctx);
        fmpz_mod_mpoly_q_clear(B, ctx);
        fmpz_mod_mpoly_q_clear(C, ctx);
        fmpz_clear(c);
        fmpz_clear(m);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
