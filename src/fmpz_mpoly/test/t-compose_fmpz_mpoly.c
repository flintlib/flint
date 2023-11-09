/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_compose_fmpz_mpoly, state)
{
    slong i, j, v;

    {
        fmpz_mpoly_t A, A1, A2, B;
        fmpz_mpoly_struct * Cp[3];
        fmpz_mpoly_struct C[3];
        fmpz_mpoly_ctx_t ctxAC, ctxB;

        fmpz_mpoly_ctx_init(ctxB, 3, ORD_LEX);
        fmpz_mpoly_ctx_init(ctxAC, 2, ORD_LEX);

        fmpz_mpoly_init(B, ctxB);
        fmpz_mpoly_init(A, ctxAC);
        fmpz_mpoly_init(A1, ctxAC);
        fmpz_mpoly_init(A2, ctxAC);
        for (i = 0; i < 3; i++)
        {
            Cp[i] = C + i;
            fmpz_mpoly_init(C + i, ctxAC);
        }

        fmpz_mpoly_set_str_pretty(B,
                "1 + x1*x2^2 + x2^9999999999999999999999999*x3^9", NULL, ctxB);

        fmpz_mpoly_set_str_pretty(C + 0, "x1 + x2", NULL, ctxAC);
        fmpz_mpoly_set_str_pretty(C + 1, "x1 - x2", NULL, ctxAC);
        fmpz_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (fmpz_mpoly_compose_fmpz_mpoly(A, B, Cp, ctxB, ctxAC) ||
            fmpz_mpoly_compose_fmpz_mpoly_horner(A1, B, Cp, ctxB, ctxAC) ||
            fmpz_mpoly_compose_fmpz_mpoly_geobucket(A2, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 1\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_set_str_pretty(C + 0, "x1", NULL, ctxAC);
        fmpz_mpoly_set_str_pretty(C + 1, "2*x2", NULL, ctxAC);
        fmpz_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (fmpz_mpoly_compose_fmpz_mpoly(A, B, Cp, ctxB, ctxAC) ||
            fmpz_mpoly_compose_fmpz_mpoly_horner(A1, B, Cp, ctxB, ctxAC) ||
            fmpz_mpoly_compose_fmpz_mpoly_geobucket(A2, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 2\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_set_str_pretty(C + 0, "2*x1", NULL, ctxAC);
        fmpz_mpoly_set_str_pretty(C + 1, "x2", NULL, ctxAC);
        fmpz_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (!fmpz_mpoly_compose_fmpz_mpoly(A, B, Cp, ctxB, ctxAC) ||
            !fmpz_mpoly_compose_fmpz_mpoly_horner(A1, B, Cp, ctxB, ctxAC) ||
            !fmpz_mpoly_compose_fmpz_mpoly_geobucket(A2, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check example 3\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(B, ctxB);
        fmpz_mpoly_clear(A, ctxAC);
        fmpz_mpoly_clear(A1, ctxAC);
        fmpz_mpoly_clear(A2, ctxAC);
        for (i = 0; i < 3; i++)
            fmpz_mpoly_clear(C + i, ctxAC);

        fmpz_mpoly_ctx_clear(ctxB);
        fmpz_mpoly_ctx_clear(ctxAC);
    }

    /* check composition with generators */
    for (i = 0; i < 20*flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctxB, ctxAC;
        fmpz_mpoly_t B, A, A1, A2, A3;
        fmpz_mpoly_struct ** C;
        slong * c;
        slong nvarsB, nvarsAC;
        slong len;
        flint_bitcnt_t exp_bits, coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctxB, state, 20);
        fmpz_mpoly_ctx_init_rand(ctxAC, state, 20);

        nvarsB = ctxB->minfo->nvars;
        nvarsAC = ctxAC->minfo->nvars;

        fmpz_mpoly_init(B, ctxB);
        fmpz_mpoly_init(A, ctxAC);
        fmpz_mpoly_init(A1, ctxAC);
        fmpz_mpoly_init(A2, ctxAC);
        fmpz_mpoly_init(A3, ctxAC);

        c = (slong *) flint_malloc(nvarsB*sizeof(slong));
        C = (fmpz_mpoly_struct **) flint_malloc(nvarsB * sizeof(fmpz_mpoly_struct *));
        for (v = 0; v < nvarsB; v++)
        {
            C[v] = (fmpz_mpoly_struct *) flint_malloc(sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(C[v], ctxAC);
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
                    fmpz_mpoly_gen(C[v], c[v], ctxAC);
                else
                    fmpz_mpoly_zero(C[v], ctxAC);
            }

            fmpz_mpoly_randtest_bits(B, state, len, coeff_bits, exp_bits, ctxB);

            fmpz_mpoly_compose_fmpz_mpoly_gen(A, B, c, ctxB, ctxAC);

            if (!fmpz_mpoly_compose_fmpz_mpoly(A1, B, C, ctxB, ctxAC) ||
                !fmpz_mpoly_compose_fmpz_mpoly_horner(A2, B, C, ctxB, ctxAC) ||
                !fmpz_mpoly_compose_fmpz_mpoly_geobucket(A3, B, C, ctxB, ctxAC))
            {
                printf("FAIL\n");
                flint_printf("Check composition success with generators\n"
                                                     "i: %wd, j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            if (!fmpz_mpoly_equal(A, A1, ctxAC) ||
                !fmpz_mpoly_equal(A, A2, ctxAC) ||
                !fmpz_mpoly_equal(A, A3, ctxAC))
            {
                printf("FAIL\n");
                flint_printf("Check composition with generators\n"
                                                     "i: %wd, j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mpoly_assert_canonical(A, ctxAC);
            fmpz_mpoly_assert_canonical(A1, ctxAC);
            fmpz_mpoly_assert_canonical(A2, ctxAC);
            fmpz_mpoly_assert_canonical(A3, ctxAC);
        }

        for (v = 0; v < nvarsB; v++)
        {
            fmpz_mpoly_clear(C[v], ctxAC);
            flint_free(C[v]);
        }
        flint_free(C);
        flint_free(c);

        fmpz_mpoly_clear(B, ctxB);
        fmpz_mpoly_clear(A, ctxAC);
        fmpz_mpoly_clear(A1, ctxAC);
        fmpz_mpoly_clear(A2, ctxAC);
        fmpz_mpoly_clear(A3, ctxAC);

        fmpz_mpoly_ctx_clear(ctxB);
        fmpz_mpoly_ctx_clear(ctxAC);
    }

    /* Check composition with identity */
    for (i = 0; i < 20*flint_test_multiplier(); i++)
    {
        slong len1;
        flint_bitcnt_t exp_bits, coeff_bits;
        fmpz_mpoly_struct ** vals1;
        fmpz_mpoly_t f, g, g1, g2;
        fmpz_mpoly_ctx_t ctx;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        vals1 = (fmpz_mpoly_struct **) flint_malloc(ctx->minfo->nvars*
                                                  sizeof(fmpz_mpoly_struct *));
        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            vals1[v] = (fmpz_mpoly_struct *) flint_malloc(sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(vals1[v], ctx);
            fmpz_mpoly_gen(vals1[v], v, ctx);
        }

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(g1, ctx);
        fmpz_mpoly_init(g2, ctx);

        len1 = n_randint(state, 200);
        exp_bits = n_randint(state, 100) + 1;
        coeff_bits = n_randint(state, 100) + 1;
        fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx);

        if (!fmpz_mpoly_compose_fmpz_mpoly(g, f, vals1, ctx, ctx) ||
            !fmpz_mpoly_compose_fmpz_mpoly_horner(g1, f, vals1, ctx, ctx) ||
            !fmpz_mpoly_compose_fmpz_mpoly_geobucket(g2, f, vals1, ctx, ctx) ||
            !fmpz_mpoly_equal(g, g1, ctx) ||
            !fmpz_mpoly_equal(g, g2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition success\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_assert_canonical(g, ctx);
        fmpz_mpoly_assert_canonical(g1, ctx);
        fmpz_mpoly_assert_canonical(g2, ctx);

        if (!fmpz_mpoly_equal(f, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition with identity\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(g1, ctx);
        fmpz_mpoly_clear(g2, ctx);

        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            fmpz_mpoly_clear(vals1[v], ctx);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check composition and evalall commute */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx1, ctx2;
        fmpz_mpoly_t f, g, g1, g2;
        fmpz_mpoly_struct ** vals1;
        fmpz_t fe, ge;
        fmpz ** vals2, ** vals3;
        slong nvars1, nvars2;
        slong len1, len2;
        slong exp_bound1, exp_bound2;
        slong coeff_bits1, coeff_bits2;

        fmpz_mpoly_ctx_init_rand(ctx1, state, 6);
        fmpz_mpoly_ctx_init_rand(ctx2, state, 6);

        nvars1 = ctx1->minfo->nvars;
        nvars2 = ctx2->minfo->nvars;

        fmpz_mpoly_init(f, ctx1);
        fmpz_mpoly_init(g, ctx2);
        fmpz_mpoly_init(g1, ctx2);
        fmpz_mpoly_init(g2, ctx2);
        fmpz_init(fe);
        fmpz_init(ge);

        len1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars1) + 1) + 1;
        len2 = n_randint(state, 10/FLINT_MAX(WORD(1), nvars2) + 1) + 1;
        exp_bound1 = n_randint(state, 15/FLINT_MAX(WORD(1), nvars1) + 2) + 2;
        exp_bound2 = n_randint(state, 15/FLINT_MAX(WORD(1), nvars2) + 2) + 3;
        coeff_bits1 = n_randint(state, 200);
        coeff_bits2 = n_randint(state, 8);

        vals1 = (fmpz_mpoly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpz_mpoly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpz_mpoly_struct *) flint_malloc(
                                                    sizeof(fmpz_mpoly_struct));
            fmpz_mpoly_init(vals1[v], ctx2);
            fmpz_mpoly_randtest_bound(vals1[v], state, len2,
                                                coeff_bits2, exp_bound2, ctx2);
        }

        vals2 = (fmpz **) flint_malloc(nvars2*sizeof(fmpz*));
        for (v = 0; v < nvars2; v++)
        {
            vals2[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals2[v]);
            fmpz_randbits(vals2[v], state, 4);
        }

        vals3 = (fmpz **) flint_malloc(nvars1*sizeof(fmpz*));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals3[v]);
            if (!fmpz_mpoly_evaluate_all_fmpz(vals3[v], vals1[v], vals2, ctx2))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits1, exp_bound1, ctx1);

        if (!fmpz_mpoly_compose_fmpz_mpoly(g, f, vals1, ctx1, ctx2) ||
            !fmpz_mpoly_compose_fmpz_mpoly_horner(g1, f, vals1, ctx1, ctx2) ||
            !fmpz_mpoly_compose_fmpz_mpoly_geobucket(g2, f, vals1, ctx1, ctx2) ||
            !fmpz_mpoly_equal(g, g1, ctx2) ||
            !fmpz_mpoly_equal(g, g2, ctx2))
        {
            printf("FAIL\n");
            flint_printf("Check composition success\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_assert_canonical(g, ctx2);
        fmpz_mpoly_assert_canonical(g1, ctx2);
        fmpz_mpoly_assert_canonical(g2, ctx2);

        if (!fmpz_mpoly_evaluate_all_fmpz(fe, f, vals3, ctx1) ||
            !fmpz_mpoly_evaluate_all_fmpz(ge, g, vals2, ctx2))
        {
            printf("FAIL\n");
            flint_printf("Check evaluation success\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_equal(fe, ge))
        {
            printf("FAIL\n");
            flint_printf("Check composition and evalall commute\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpz_mpoly_clear(vals1[v], ctx2);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        for (v = 0; v < nvars2; v++)
        {
            fmpz_clear(vals2[v]);
            flint_free(vals2[v]);
        }
        flint_free(vals2);

        for (v = 0; v < nvars1; v++)
        {
            fmpz_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        fmpz_mpoly_clear(f, ctx1);
        fmpz_mpoly_clear(g, ctx2);
        fmpz_mpoly_clear(g1, ctx2);
        fmpz_mpoly_clear(g2, ctx2);

        fmpz_clear(fe);
        fmpz_clear(ge);

        fmpz_mpoly_ctx_clear(ctx1);
        fmpz_mpoly_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
