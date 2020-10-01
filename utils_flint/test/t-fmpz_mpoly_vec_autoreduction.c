/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "calcium.h"
#include "utils_flint.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("fmpz_mpoly_vec_autoreduction...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * calcium_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_vec_t F, G;
        slong nvars;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_vec_init(F, ctx);
        fmpz_mpoly_vec_init(G, ctx);

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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
