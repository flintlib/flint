/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_mul, state)
{
    fmpz_mod_mat_t A, B, B1, B2, C, C1, C2, D;
    slong max_threads = 5;
    slong i;

    /* test A*(B1+B1) = A*B1 + A*B2 */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong m, n, k;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        if (n_randint(state, 100) == 0)
        {
            m = n_randint(state, 300) + 100;
            n = n_randint(state, 300) + 100;
            k = n_randint(state, 300) + 100;
        }
        else if (n_randint(state, 10) == 0)
        {
            m = n_randint(state, 50);
            n = n_randint(state, 50);
            k = n_randint(state, 50);
        }
        else
        {
            m = n_randint(state, 8);
            n = n_randint(state, 8);
            k = n_randint(state, 8);
        }

        fmpz_mod_mat_init(A, m, n, ctx);
        fmpz_mod_mat_init(B, n, k, ctx);
        fmpz_mod_mat_init(B1, n, k, ctx);
        fmpz_mod_mat_init(B2, n, k, ctx);
        fmpz_mod_mat_init(C, m, k, ctx);
        fmpz_mod_mat_init(C1, m, k, ctx);
        fmpz_mod_mat_init(C2, m, k, ctx);
        fmpz_mod_mat_init(D, m, k, ctx);

        fmpz_mod_mat_randtest(A, state, ctx);
        fmpz_mod_mat_randtest(B1, state, ctx);
        fmpz_mod_mat_randtest(B2, state, ctx);

        /* Make sure noise in the output is ok */
        fmpz_mod_mat_randtest(C, state, ctx);
        fmpz_mod_mat_randtest(C1, state, ctx);
        fmpz_mod_mat_randtest(C2, state, ctx);

        fmpz_mod_mat_mul(C1, A, B1, ctx);
        fmpz_mod_mat_mul(C2, A, B2, ctx);
        fmpz_mod_mat_add(B, B1, B2, ctx);
        fmpz_mod_mat_mul(C, A, B, ctx);
        fmpz_mod_mat_add(D, C1, C2, ctx);

        if (!fmpz_mod_mat_equal(C, D, ctx))
        {
            flint_printf("FAIL: results not equal\n\n");
            fmpz_mod_mat_print_pretty(A, ctx); flint_printf("\n\n");
            fmpz_mod_mat_print_pretty(B1, ctx); flint_printf("\n\n");
            fmpz_mod_mat_print_pretty(B2, ctx); flint_printf("\n\n");
            fmpz_mod_mat_print_pretty(C, ctx); flint_printf("\n\n");
            fmpz_mod_mat_print_pretty(D, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (n == k)
        {
            fmpz_mod_mat_mul(A, A, B, ctx);

            if (!fmpz_mod_mat_equal(A, C, ctx))
            {
                flint_printf("FAIL: aliasing failed\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(B1, ctx);
        fmpz_mod_mat_clear(B2, ctx);
        fmpz_mod_mat_clear(C, ctx);
        fmpz_mod_mat_clear(C1, ctx);
        fmpz_mod_mat_clear(C2, ctx);
        fmpz_mod_mat_clear(D, ctx);
        fmpz_mod_ctx_clear(ctx);
    }

    /* Test aliasing with windows */
    {
        fmpz_mod_mat_t A_window;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_ui(ctx, 3);

        fmpz_mod_mat_init(A, 2, 2, ctx);
        fmpz_mod_mat_init(B, 2, 2, ctx);

        fmpz_mod_mat_window_init(A_window, A, 0, 0, 2, 2, ctx);

        fmpz_mod_mat_one(A, ctx);
        fmpz_mod_mat_one(B, ctx);
        fmpz_set_ui(fmpz_mod_mat_entry(B, 0, 1), 1);
        fmpz_set_ui(fmpz_mod_mat_entry(B, 1, 0), 1);

        fmpz_mod_mat_mul(A_window, B, A_window, ctx);

        if (!fmpz_mod_mat_equal(A, B, ctx))
        {
            flint_printf("FAIL: window aliasing failed\n");
            fmpz_mod_mat_print_pretty(A, ctx); flint_printf("\n\n");
            fmpz_mod_mat_print_pretty(B, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_window_clear(A_window, ctx);
        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
