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

    flint_printf("div_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * calcium_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_q_t A, B, C;
        fmpz_t c;

        fmpz_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX);

        fmpz_mpoly_q_init(A, ctx);
        fmpz_mpoly_q_init(B, ctx);
        fmpz_mpoly_q_init(C, ctx);
        fmpz_init(c);

        fmpz_mpoly_q_randtest(A, state, 10, 2 + n_randint(state, 100), 5, ctx);
        fmpz_randtest_not_zero(c, state, 10);

        fmpz_mpoly_q_div_fmpz(B, A, c, ctx);

        fmpz_mpoly_set(fmpz_mpoly_q_numref(C), fmpz_mpoly_q_numref(A), ctx);
        fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_q_denref(C), fmpz_mpoly_q_denref(A), c, ctx);
        fmpz_mpoly_q_canonicalise(C, ctx);

        if (!fmpz_mpoly_q_equal(B, C, ctx))
        {
            flint_printf("FAIL\n");
            flint_printf("A = "); fmpz_mpoly_q_print_pretty(A, NULL, ctx); flint_printf("\n\n");
            flint_printf("B = "); fmpz_mpoly_q_print_pretty(B, NULL, ctx); flint_printf("\n\n");
            flint_printf("C = "); fmpz_mpoly_q_print_pretty(C, NULL, ctx); flint_printf("\n\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mpoly_q_clear(A, ctx);
        fmpz_mpoly_q_clear(B, ctx);
        fmpz_mpoly_q_clear(C, ctx);
        fmpz_clear(c);

        fmpz_mpoly_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

