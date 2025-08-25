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

TEST_FUNCTION_START(fmpz_mod_mpoly_q_add, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_q_t A, B, C, D;
        fmpz_mod_mpoly_t t, u;
        fmpz_t m;

        fmpz_init(m);
        fmpz_randtest_unsigned(m, state, n_randint(state, 2) ? 4 : n_randint(state, 100));
        fmpz_nextprime(m, m, 0);
        fmpz_mod_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX, m);

        fmpz_mod_mpoly_q_init(A, ctx);
        fmpz_mod_mpoly_q_init(B, ctx);
        fmpz_mod_mpoly_q_init(C, ctx);
        fmpz_mod_mpoly_q_init(D, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);

        fmpz_mod_mpoly_q_randtest(A, state, 4, 4, ctx);
        fmpz_mod_mpoly_q_randtest(B, state, 4, 4, ctx);

        fmpz_mod_mpoly_q_add(C, A, B, ctx);

        fmpz_mod_mpoly_mul(t, fmpz_mod_mpoly_q_numref(A), fmpz_mod_mpoly_q_denref(B), ctx);
        fmpz_mod_mpoly_mul(u, fmpz_mod_mpoly_q_numref(B), fmpz_mod_mpoly_q_denref(A), ctx);
        fmpz_mod_mpoly_add(fmpz_mod_mpoly_q_numref(D), t, u, ctx);
        fmpz_mod_mpoly_mul(fmpz_mod_mpoly_q_denref(D), fmpz_mod_mpoly_q_denref(A), fmpz_mod_mpoly_q_denref(B), ctx);

        fmpz_mod_mpoly_q_canonicalise(D, ctx);

        if (!fmpz_mod_mpoly_q_equal(C, D, ctx))
        {
            flint_printf("FAIL\n");
            flint_printf("A = "); fmpz_mod_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mod_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mod_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("D = "); fmpz_mod_mpoly_q_print_pretty(D, NULL, ctx); flint_printf("\n\n");
            flint_printf("mod = ");fmpz_print(ctx->ffinfo->n);flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            fmpz_mod_mpoly_q_set(C, A, ctx);
            fmpz_mod_mpoly_q_add(C, C, B, ctx);
        }
        else
        {
            fmpz_mod_mpoly_q_set(C, B, ctx);
            fmpz_mod_mpoly_q_add(C, A, C, ctx);
        }

        if (!fmpz_mod_mpoly_q_equal(C, D, ctx))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("A = "); fmpz_mod_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mod_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mod_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("D = "); fmpz_mod_mpoly_q_print_pretty(D, NULL, ctx); flint_printf("\n\n");
            flint_printf("mod = ");fmpz_print(ctx->ffinfo->n);flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mod_mpoly_q_clear(A, ctx);
        fmpz_mod_mpoly_q_clear(B, ctx);
        fmpz_mod_mpoly_q_clear(C, ctx);
        fmpz_mod_mpoly_q_clear(D, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(u, ctx);
        fmpz_clear(m);

        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
