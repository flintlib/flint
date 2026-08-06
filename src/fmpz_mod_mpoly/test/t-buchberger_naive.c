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
#include "calcium.h"
#include "fmpz_mod_mpoly.h"
#include "mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_buchberger_naive, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_vec_t F, G, H;
        slong nvars;
        fmpz_t m;

        fmpz_init(m);
        fmpz_randtest_unsigned(m, state, n_randint(state, 2) ? 4 : n_randint(state, 100));
        fmpz_nextprime(m, m, 0);
        fmpz_mod_mpoly_ctx_init(ctx, 1 + n_randint(state, 4), ORD_LEX, m);
        nvars = ctx->minfo->nvars;

        fmpz_mod_mpoly_vec_init(F, 0, ctx);
        fmpz_mod_mpoly_vec_init(G, 0, ctx);
        fmpz_mod_mpoly_vec_init(H, 0, ctx);

        
        // flint_printf("iter %wd   %wd  %d    %{fmpz}\n\n", iter, nvars, ctx->minfo->ord, m);
        // printf("--------------------------------------------------------------------\n");
        

        if (nvars == 4)
            fmpz_mod_mpoly_vec_randtest_not_zero(F, state, 1 + n_randint(state, 3), 1 + n_randint(state, 3), 1 + n_randint(state, 2), ctx);
        else if (nvars == 3)
            fmpz_mod_mpoly_vec_randtest_not_zero(F, state, 1 + n_randint(state, 4), 1 + n_randint(state, 4), 1 + n_randint(state, 2), ctx);
        else
            fmpz_mod_mpoly_vec_randtest_not_zero(F, state, 1 + n_randint(state, 5), 10+fmpz_bits(m), 1 + n_randint(state, 3), ctx);

        // flint_printf("F = "); fmpz_mod_mpoly_vec_print(F, ctx); flint_printf("\n");

        fmpz_mod_mpoly_buchberger_naive(G, F, ctx);

        // flint_printf("G = "); fmpz_mod_mpoly_vec_print(G, ctx); flint_printf("\n");

        if (!fmpz_mod_mpoly_vec_is_groebner(G, F, ctx))
        {
            flint_printf("FAIL\n\n");
            mpoly_ordering_print(ctx->minfo->ord); printf("\n");
            flint_printf("F = "); fmpz_mod_mpoly_vec_print(F, ctx); flint_printf("\n");
            flint_printf("G = "); fmpz_mod_mpoly_vec_print(G, ctx); flint_printf("\n");
            flint_abort();
        }
        
        fmpz_mod_mpoly_vec_autoreduction_groebner(H, G, ctx);

        if (!fmpz_mod_mpoly_vec_is_groebner(H, F, ctx) || !fmpz_mod_mpoly_vec_is_autoreduced(H, ctx))
        {
            flint_printf("FAIL (reduced GB)\n\n");
            mpoly_ordering_print(ctx->minfo->ord); printf("\n");
            flint_printf("F = "); fmpz_mod_mpoly_vec_print(F, ctx); flint_printf("\n");
            flint_printf("G = "); fmpz_mod_mpoly_vec_print(G, ctx); flint_printf("\n");
            flint_printf("H = "); fmpz_mod_mpoly_vec_print(H, ctx); flint_printf("\n");
            flint_abort();
        }

        fmpz_mod_mpoly_vec_clear(F, ctx);
        fmpz_mod_mpoly_vec_clear(G, ctx);
        fmpz_mod_mpoly_vec_clear(H, ctx);

        fmpz_mod_mpoly_ctx_clear(ctx);

        fmpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
