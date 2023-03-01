/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("add....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * calcium_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_t A, B, C, D;
        fmpz_mpoly_t t, u;

        fmpz_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX);

        fmpz_mpoly_q_init(A, ctx);
        fmpz_mpoly_q_init(B, ctx);
        fmpz_mpoly_q_init(C, ctx);
        fmpz_mpoly_q_init(D, ctx);
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_init(u, ctx);

        fmpz_mpoly_q_randtest(A, state, 5, 2 + n_randint(state, 100), 5, ctx);
        fmpz_mpoly_q_randtest(B, state, 5, 2 + n_randint(state, 100), 5, ctx);

        fmpz_mpoly_q_add(C, A, B, ctx);

        fmpz_mpoly_mul(t, fmpz_mpoly_q_numref(A), fmpz_mpoly_q_denref(B), ctx);
        fmpz_mpoly_mul(u, fmpz_mpoly_q_numref(B), fmpz_mpoly_q_denref(A), ctx);
        fmpz_mpoly_add(fmpz_mpoly_q_numref(D), t, u, ctx);
        fmpz_mpoly_mul(fmpz_mpoly_q_denref(D), fmpz_mpoly_q_denref(A), fmpz_mpoly_q_denref(B), ctx);

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
            fmpz_mpoly_q_add(C, C, B, ctx);
        }
        else
        {
            fmpz_mpoly_q_set(C, B, ctx);
            fmpz_mpoly_q_add(C, A, C, ctx);
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
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(u, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

