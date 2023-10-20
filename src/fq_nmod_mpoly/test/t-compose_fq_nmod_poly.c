/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_compose_fq_nmod_poly, state)
{
    slong i, v;

    {
        fq_nmod_poly_t A;
        fq_nmod_mpoly_t B;
        fq_nmod_poly_struct * Cp[3];
        fq_nmod_poly_struct C[3];
        fq_nmod_mpoly_ctx_t ctxB;
        fmpz one = 1, two = 2, three = 3;

        fq_nmod_mpoly_ctx_init_deg(ctxB, 3, ORD_LEX, 13, 3);

        fq_nmod_mpoly_init(B, ctxB);
        fq_nmod_poly_init(A, ctxB->fqctx);
        for (i = 0; i < 3; i++)
        {
            Cp[i] = C + i;
            fq_nmod_poly_init(C + i, ctxB->fqctx);
        }

        fq_nmod_mpoly_set_str_pretty(B,
                "1 + x1*x2^2 + x2^9999999999999999999999999*x3^9", NULL, ctxB);

        fq_nmod_poly_zero(C + 0, ctxB->fqctx);
        fq_nmod_poly_zero(C + 1, ctxB->fqctx);
        fq_nmod_poly_zero(C + 2, ctxB->fqctx);
        fq_nmod_poly_set_coeff_fmpz(C + 0, 1, &one, ctxB->fqctx);
        fq_nmod_poly_set_coeff_fmpz(C + 1, 2, &two, ctxB->fqctx);
        fq_nmod_poly_set_coeff_fmpz(C + 2, 3, &three, ctxB->fqctx);
        if (fq_nmod_mpoly_compose_fq_nmod_poly(A, B, Cp, ctxB))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 1\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_poly_zero(C + 0, ctxB->fqctx);
        fq_nmod_poly_zero(C + 1, ctxB->fqctx);
        fq_nmod_poly_zero(C + 2, ctxB->fqctx);
        fq_nmod_poly_set_coeff_fmpz(C + 0, 0, &one, ctxB->fqctx);
        fq_nmod_poly_set_coeff_fmpz(C + 1, 0, &two, ctxB->fqctx);
        fq_nmod_poly_set_coeff_fmpz(C + 2, 0, &three, ctxB->fqctx);
        if (!fq_nmod_mpoly_compose_fq_nmod_poly(A, B, Cp, ctxB))
        {
            printf("FAIL\n");
            flint_printf("Check example 2\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (!fq_nmod_poly_is_zero(A, ctxB->fqctx))
        {
            printf("FAIL\n");
            flint_printf("Check example 2 equality\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(B, ctxB);
        fq_nmod_poly_clear(A, ctxB->fqctx);
        for (i = 0; i < 3; i++)
            fq_nmod_poly_clear(C + i, ctxB->fqctx);

        fq_nmod_mpoly_ctx_clear(ctxB);
    }

    /* Check composition and evalall commute */
    for (i = 0; i < 20*flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f;
        fq_nmod_poly_t g;
        fq_nmod_poly_struct ** vals1;
        fq_nmod_t fe, ge;
        fq_nmod_t vals2;
        fq_nmod_struct ** vals3;
        slong nvars1;
        slong len1, len2;
        slong exp_bound1;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        nvars1 = ctx->minfo->nvars;

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        fq_nmod_init(fe, ctx->fqctx);
        fq_nmod_init(ge, ctx->fqctx);

        len1 = n_randint(state, 40/FLINT_MAX(WORD(1), nvars1) + 1);
        len2 = n_randint(state, 20);
        exp_bound1 = n_randint(state, 100/FLINT_MAX(WORD(1), nvars1) + 2) + 1;

        fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);

        vals1 = (fq_nmod_poly_struct **) flint_malloc(nvars1*
                                                sizeof(fq_nmod_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fq_nmod_poly_struct *) flint_malloc(
                                                  sizeof(fq_nmod_poly_struct));
            fq_nmod_poly_init(vals1[v], ctx->fqctx);
            fq_nmod_poly_randtest(vals1[v], state, len2, ctx->fqctx);
        }

        fq_nmod_init(vals2, ctx->fqctx);
        fq_nmod_randtest(vals2, state, ctx->fqctx);

        vals3 = (fq_nmod_struct **) flint_malloc(nvars1*sizeof(fq_nmod_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fq_nmod_struct *) flint_malloc(sizeof(fq_nmod_struct));
            fq_nmod_init(vals3[v], ctx->fqctx);
            fq_nmod_poly_evaluate_fq_nmod(vals3[v], vals1[v], vals2, ctx->fqctx);
        }

        if (fq_nmod_mpoly_total_degree_si(f, ctx) < 50)
        {
            if (!fq_nmod_mpoly_compose_fq_nmod_poly(g, f, vals1, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check composition success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_evaluate_all_fq_nmod(fe, f, vals3, ctx);
            fq_nmod_poly_evaluate_fq_nmod(ge, g, vals2, ctx->fqctx);

            if (!fq_nmod_equal(fe, ge, ctx->fqctx))
            {
                printf("FAIL\n");
                flint_printf("Check composition and evalall commute\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars1; v++)
        {
            fq_nmod_poly_clear(vals1[v], ctx->fqctx);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fq_nmod_clear(vals2, ctx->fqctx);

        for (v = 0; v < nvars1; v++)
        {
            fq_nmod_clear(vals3[v], ctx->fqctx);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_poly_clear(g, ctx->fqctx);
        fq_nmod_clear(fe, ctx->fqctx);
        fq_nmod_clear(ge, ctx->fqctx);

        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
