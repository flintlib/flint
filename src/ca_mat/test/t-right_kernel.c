/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("right_kernel...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, X, AX;
        slong r, c, rank, nullity;
        int success, success2;

        ca_ctx_init(ctx);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        ca_mat_init(A, r, c, ctx);
        ca_mat_init(X, r, c, ctx);

        ca_mat_randtest(A, state, 1, 5, ctx);

        success = ca_mat_right_kernel(X, A, ctx);

        if (success)
        {
            success2 = ca_mat_rank(&rank, A, ctx);

            if (success2)
            {
                nullity = c - rank;

                if (nullity != ca_mat_ncols(X))
                {
                    flint_printf("FAIL:\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); ca_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank, nullity);
                    flint_abort();
                }

                ca_mat_init(AX, r, nullity, ctx);

                ca_mat_mul(AX, A, X, ctx);

                if (ca_mat_check_is_zero(AX, ctx) == T_FALSE)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); ca_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("AX: "); ca_mat_print(AX, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank, nullity);
                    flint_abort();
                }

                success = ca_mat_rank(&nullity, X, ctx);

                if (!success || (nullity != ca_mat_ncols(X)))
                {
                    flint_printf("FAIL (2):\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("X: "); ca_mat_print(X, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank, nullity);
                    flint_abort();
                }

                ca_mat_clear(AX, ctx);
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(X, ctx);

        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
