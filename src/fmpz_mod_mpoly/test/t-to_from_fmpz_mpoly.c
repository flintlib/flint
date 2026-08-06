/*
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
#include "fmpz_mpoly.h"
#include "mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_get_fmpz_mpoly, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_ctx_t ctxm;
        fmpz_mod_mpoly_t A, B;
        fmpz_mpoly_t C, D;
        slong nvars, i;
        fmpz_t m, c;

        nvars = 1 + n_randint(state, 3);
        fmpz_init(m);
        fmpz_init(c);
        fmpz_randtest_unsigned(m, state, n_randint(state, 2) ? 4 : n_randint(state, 90));
        fmpz_nextprime(m, m, 0);
        fmpz_mod_mpoly_ctx_init(ctxm, nvars, ORD_LEX, m);
        fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);

        fmpz_mod_mpoly_init(A, ctxm);
        fmpz_mod_mpoly_init(B, ctxm);
        fmpz_mpoly_init(C, ctx);
        fmpz_mpoly_init(D, ctx);

        fmpz ** exp = (fmpz **) flint_malloc(ctx->minfo->nvars*sizeof(fmpz *));

        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            exp[i] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(exp[i]);
        }

        // flint_printf("iter %wd   %wd  %d    %{fmpz}\n\n", iter, nvars, ctx->minfo->ord, m);
        // printf("--------------------------------------------------------------------\n");
        

        fmpz_mpoly_randtest_bound(C, state, 1 + n_randint(state, 4), 1 + n_randint(state, 4), 1 + n_randint(state, 4), ctx);

        fmpz_mod_mpoly_resize(B, C->length, ctxm);

        for (i = 0; i< C->length; i++)
        {
            fmpz_mpoly_get_term_coeff_fmpz(c, C, i, ctx);
            fmpz_mpoly_get_term_exp_fmpz(exp, C, i, ctx);
            
            fmpz_mod_mpoly_set_term_exp_fmpz(B, i, exp, ctxm);
            fmpz_mod_mpoly_set_term_coeff_fmpz(B, i, c, ctxm);
        }

        fmpz_mod_mpoly_combine_like_terms(B, ctxm);
        fmpz_mod_mpoly_set_fmpz_mpoly(A, C, ctxm, ctx);

        if (!fmpz_mod_mpoly_equal(A, B, ctxm))
        {
            flint_printf("FAIL\n\n");
            mpoly_ordering_print(ctx->minfo->ord); printf("\n");
            flint_printf("C = "); fmpz_mpoly_print_pretty(C, NULL, ctx); flint_printf("\n");
            flint_printf("A = "); fmpz_mod_mpoly_print_pretty(A, NULL, ctxm); flint_printf("\n");
            flint_printf("B = "); fmpz_mod_mpoly_print_pretty(B, NULL, ctxm); flint_printf("\n");
            flint_abort();
        }

        fmpz_mpoly_resize(C, A->length, ctx);

        for (i = 0; i< A->length; i++)
        {
            fmpz_mod_mpoly_get_term_coeff_fmpz(c, A, i, ctxm);
            fmpz_mod_mpoly_get_term_exp_fmpz(exp, A, i, ctxm);
            
            fmpz_mpoly_set_term_exp_fmpz(C, i, exp, ctx);
            fmpz_mpoly_set_term_coeff_fmpz(C, i, c, ctx);
        }

        fmpz_mod_mpoly_get_fmpz_mpoly(D, A, ctx);  

        if (!fmpz_mpoly_equal(C, D, ctx))
        {
            flint_printf("FAIL\n\n");
            mpoly_ordering_print(ctx->minfo->ord); printf("\n");
            flint_printf("C = "); fmpz_mpoly_print_pretty(C, NULL, ctx); flint_printf("\n");
            flint_printf("D = "); fmpz_mpoly_print_pretty(D, NULL, ctx); flint_printf("\n");
            flint_printf("A = "); fmpz_mod_mpoly_print_pretty(A, NULL, ctxm); flint_printf("\n");
            flint_abort();
        }

        for (i = 0; i < ctx->minfo->nvars; i++)
        {
            fmpz_clear(exp[i]);
            flint_free(exp[i]);
        }
        flint_free(exp);

        fmpz_mod_mpoly_clear(A, ctxm);
        fmpz_mod_mpoly_clear(B, ctxm);
        fmpz_mpoly_clear(C, ctx);
        fmpz_mpoly_clear(D, ctx);

        fmpz_mod_mpoly_ctx_clear(ctxm);
        fmpz_mpoly_ctx_clear(ctx);

        fmpz_clear(m);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
