/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mat.h"
#include "fq.h"
#include "fq_embed.h"

TEST_FUNCTION_START(fq_embed_mono_dual_matrix, state)
{
    int i;

    /* Check that the two functions are inverse of one another */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fmpz_mod_mat_t m2d, d2m, one, two;
        slong d;

        fq_ctx_init_randtest(ctx, state, 0);
        d = fq_ctx_degree(ctx);

        fmpz_mod_mat_init(m2d, d, d, ctx->ctxp);
        fmpz_mod_mat_init(d2m, d, d, ctx->ctxp);
        fmpz_mod_mat_init(one, d, d, ctx->ctxp);
        fmpz_mod_mat_init(two, d, d, ctx->ctxp);

        fq_embed_mono_to_dual_matrix(m2d, ctx);
        fq_embed_dual_to_mono_matrix(d2m, ctx);
        fmpz_mod_mat_mul(one, m2d, d2m, ctx->ctxp);

        fmpz_mod_mat_one(two, ctx->ctxp);

        if (!fmpz_mod_mat_equal(one, two, ctx->ctxp))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_ctx_print(ctx), flint_printf("\n");
            flint_printf("Mono -> Dual\n"),
                fmpz_mod_mat_print_pretty(m2d, ctx->ctxp), flint_printf("\nDual -> Mono\n"),
                fmpz_mod_mat_print_pretty(d2m, ctx->ctxp), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(m2d, ctx->ctxp);
        fmpz_mod_mat_clear(d2m, ctx->ctxp);
        fmpz_mod_mat_clear(one, ctx->ctxp);
        fmpz_mod_mat_clear(two, ctx->ctxp);
        fq_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
