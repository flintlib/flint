/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq_mpoly.h"

int
main(void)
{
    slong i, j, v;
    FLINT_TEST_INIT(state);

    flint_printf("compose_fmpq_mpoly....");
    fflush(stdout);

    {
        fmpq_mpoly_t A, A1, A2, B;
        fmpq_mpoly_struct * Cp[3];
        fmpq_mpoly_struct C[3];
        fmpq_mpoly_ctx_t ctxAC, ctxB;

        fmpq_mpoly_ctx_init(ctxB, 3, ORD_LEX);
        fmpq_mpoly_ctx_init(ctxAC, 2, ORD_LEX);

        fmpq_mpoly_init(B, ctxB);
        fmpq_mpoly_init(A, ctxAC);
        fmpq_mpoly_init(A1, ctxAC);
        fmpq_mpoly_init(A2, ctxAC);
        for (i = 0; i < 3; i++)
        {
            Cp[i] = C + i;
            fmpq_mpoly_init(C + i, ctxAC);
        }

        fmpq_mpoly_set_str_pretty(B,
                "1 + x1*x2^2 + x2^9999999999999999999999999*x3^9", NULL, ctxB);

        fmpq_mpoly_set_str_pretty(C + 0, "x1 + x2", NULL, ctxAC);
        fmpq_mpoly_set_str_pretty(C + 1, "x1 - x2", NULL, ctxAC);
        fmpq_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (fmpq_mpoly_compose_fmpq_mpoly(A, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 1\n", i);
            flint_abort();
        }

        fmpq_mpoly_set_str_pretty(C + 0, "x1", NULL, ctxAC);
        fmpq_mpoly_set_str_pretty(C + 1, "2*x2", NULL, ctxAC);
        fmpq_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (fmpq_mpoly_compose_fmpq_mpoly(A, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 2\n", i);
            flint_abort();
        }

        fmpq_mpoly_set_str_pretty(C + 0, "2*x1", NULL, ctxAC);
        fmpq_mpoly_set_str_pretty(C + 1, "x2", NULL, ctxAC);
        fmpq_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (!fmpq_mpoly_compose_fmpq_mpoly(A, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check example 3\n", i);
            flint_abort();
        }

        fmpq_mpoly_clear(B, ctxB);
        fmpq_mpoly_clear(A, ctxAC);
        fmpq_mpoly_clear(A1, ctxAC);
        fmpq_mpoly_clear(A2, ctxAC);
        for (i = 0; i < 3; i++)
            fmpq_mpoly_clear(C + i, ctxAC);

        fmpq_mpoly_ctx_clear(ctxB);
        fmpq_mpoly_ctx_clear(ctxAC);
    }

    /* check composition with generators */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctxB, ctxAC;
        fmpq_mpoly_t B, A, A1;
        fmpq_mpoly_struct ** C;
        slong * c;
        slong nvarsB, nvarsAC;
        slong len;
        flint_bitcnt_t exp_bits, coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctxB, state, 20);
        fmpq_mpoly_ctx_init_rand(ctxAC, state, 20);

        nvarsB = ctxB->zctx->minfo->nvars;
        nvarsAC = ctxAC->zctx->minfo->nvars;

        fmpq_mpoly_init(B, ctxB);
        fmpq_mpoly_init(A, ctxAC);
        fmpq_mpoly_init(A1, ctxAC);

        c = (slong *) flint_malloc(nvarsB*sizeof(slong));
        C = (fmpq_mpoly_struct **) flint_malloc(nvarsB * sizeof(fmpq_mpoly_struct *));
        for (v = 0; v < nvarsB; v++)
        {
            C[v] = (fmpq_mpoly_struct *) flint_malloc(sizeof(fmpq_mpoly_struct));
            fmpq_mpoly_init(C[v], ctxAC);
        }

        for (j = 0; j < 4; j++)
        {
            len = n_randint(state, 200);
            exp_bits = n_randint(state, 100) + 1;
            coeff_bits = n_randint(state, 100) + 1;

            for (v = 0; v < nvarsB; v++)
            {
                c[v] = n_randint(state, nvarsAC + 2) - 2;
                if (c[v] >= 0)
                    fmpq_mpoly_gen(C[v], c[v], ctxAC);
                else
                    fmpq_mpoly_zero(C[v], ctxAC);
            }

            fmpq_mpoly_randtest_bits(B, state, len, coeff_bits, exp_bits, ctxB);

            fmpq_mpoly_compose_fmpq_mpoly_gen(A, B, c, ctxB, ctxAC);

            if (!fmpq_mpoly_compose_fmpq_mpoly(A1, B, C, ctxB, ctxAC))
            {
                printf("FAIL\n");
                flint_printf("Check composition success with generators\n"
                                                     "i: %wd, j: %wd\n", i, j);
                flint_abort();
            }

            if (!fmpq_mpoly_equal(A, A1, ctxAC))
            {
                printf("FAIL\n");
                flint_printf("Check composition with generators\n"
                                                     "i: %wd, j: %wd\n", i, j);
                flint_abort();
            }

            fmpq_mpoly_assert_canonical(A, ctxAC);
            fmpq_mpoly_assert_canonical(A1, ctxAC);
        }

        for (v = 0; v < nvarsB; v++)
        {
            fmpq_mpoly_clear(C[v], ctxAC);
            flint_free(C[v]);
        }
        flint_free(C);
        flint_free(c);

        fmpq_mpoly_clear(B, ctxB);
        fmpq_mpoly_clear(A, ctxAC);
        fmpq_mpoly_clear(A1, ctxAC);

        fmpq_mpoly_ctx_clear(ctxB);
        fmpq_mpoly_ctx_clear(ctxAC);
    }

    /* Check composition with identity */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        slong nvars, len1;
        flint_bitcnt_t exp_bits, coeff_bits;
        fmpq_mpoly_struct ** vals1;
        fmpq_mpoly_t f, g;
        fmpq_mpoly_ctx_t ctx;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->zctx->minfo->nvars;

        vals1 = (fmpq_mpoly_struct **) flint_malloc(nvars*sizeof(fmpq_mpoly_struct *));
        for (v = 0; v < nvars; v++)
        {
            vals1[v] = (fmpq_mpoly_struct *) flint_malloc(sizeof(fmpq_mpoly_struct)); 
            fmpq_mpoly_init(vals1[v], ctx);
            fmpq_mpoly_gen(vals1[v], v, ctx);            
        }

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);

        len1 = n_randint(state, 200);
        exp_bits = n_randint(state, 100) + 1;
        coeff_bits = n_randint(state, 100) + 1;
        fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx);
        if (!fmpq_mpoly_compose_fmpq_mpoly(g, f, vals1, ctx, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition success\ni: %wd\n", i);
            flint_abort();
        }

        if (!fmpq_mpoly_equal(f, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition with identity\ni: %wd\n", i);
            flint_abort();
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(f, ctx);

        for (v = 0; v < nvars; v++)
        {
            fmpq_mpoly_clear(vals1[v], ctx);
            flint_free(vals1[v]);            
        }
        flint_free(vals1);
    }

    /* Check composition and evalall commute */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx1, ctx2;
        fmpq_mpoly_t f, g;
        fmpq_mpoly_struct ** vals1;
        fmpq_t fe, ge;
        fmpq ** vals2, ** vals3;
        slong nvars1, nvars2;
        slong len1, len2;
        slong exp_bound1, exp_bound2;
        slong coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx1, state, 10);
        fmpq_mpoly_ctx_init_rand(ctx2, state, 10);
        nvars1 = ctx1->zctx->minfo->nvars;
        nvars2 = ctx2->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx1);
        fmpq_mpoly_init(g, ctx2);
        fmpq_init(fe);
        fmpq_init(ge);

        len1 = n_randint(state, 70/nvars1 + 1);
        len2 = n_randint(state, 25/nvars2 + 1);
        exp_bound1 = n_randint(state, 15/nvars1 + 2) + 1;
        exp_bound2 = n_randint(state, 15/nvars2 + 2) + 1;
        coeff_bits = n_randint(state, 10);

        vals1 = (fmpq_mpoly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpq_mpoly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpq_mpoly_struct *) flint_malloc(
                                                    sizeof(fmpq_mpoly_struct)); 
            fmpq_mpoly_init(vals1[v], ctx2);
            fmpq_mpoly_randtest_bound(vals1[v], state, len2,
                                                 coeff_bits, exp_bound2, ctx2);
        }

        vals2 = (fmpq **) flint_malloc(nvars2*sizeof(fmpq*));
        for (v = 0; v < nvars2; v++)
        {
            vals2[v] = (fmpq *) flint_malloc(sizeof(fmpq));
            fmpq_init(vals2[v]);
            fmpq_randbits(vals2[v], state, 4);
        }

        vals3 = (fmpq **) flint_malloc(nvars1*sizeof(fmpq*));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpq *) flint_malloc(sizeof(fmpq)); 
            fmpq_init(vals3[v]);
            if (!fmpq_mpoly_evaluate_all_fmpq(vals3[v], vals1[v], vals2, ctx2))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                flint_abort();
            }
        }

        fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx1);

        if (!fmpq_mpoly_compose_fmpq_mpoly(g, f, vals1, ctx1, ctx2))
        {
            printf("FAIL\n");
            flint_printf("Check composition success\ni: %wd\n", i);
            flint_abort();
        }
        fmpq_mpoly_assert_canonical(g, ctx2);

        if (!fmpq_mpoly_evaluate_all_fmpq(fe, f, vals3, ctx1) ||
            !fmpq_mpoly_evaluate_all_fmpq(ge, g, vals2, ctx2))
        {
            printf("FAIL\n");
            flint_printf("Check evaluation success\ni: %wd\n", i);
            flint_abort();
        }

        if (!fmpq_equal(fe, ge))
        {
            printf("FAIL\n");
            flint_printf("Check composition and evalall commute\ni: %wd\n", i);
            flint_abort();
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpq_mpoly_clear(vals1[v], ctx2);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        for (v = 0; v < nvars2; v++)
        {
            fmpq_clear(vals2[v]);
            flint_free(vals2[v]);
        }
        flint_free(vals2);

        for (v = 0; v < nvars1; v++)
        {
            fmpq_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        fmpq_mpoly_clear(f, ctx1);
        fmpq_mpoly_clear(g, ctx2);

        fmpq_clear(fe);
        fmpq_clear(ge);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
