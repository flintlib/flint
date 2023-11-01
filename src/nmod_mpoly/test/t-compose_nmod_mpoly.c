/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_compose_nmod_mpoly, state)
{
    slong i, j, v;

    {
        nmod_mpoly_t A, A1, A2, B;
        nmod_mpoly_struct * Cp[3];
        nmod_mpoly_struct C[3];
        nmod_mpoly_ctx_t ctxAC, ctxB;

        nmod_mpoly_ctx_init(ctxB, 3, ORD_LEX, 13);
        nmod_mpoly_ctx_init(ctxAC, 2, ORD_LEX, 13);

        nmod_mpoly_init(B, ctxB);
        nmod_mpoly_init(A, ctxAC);
        nmod_mpoly_init(A1, ctxAC);
        nmod_mpoly_init(A2, ctxAC);
        for (i = 0; i < 3; i++)
        {
            Cp[i] = C + i;
            nmod_mpoly_init(C + i, ctxAC);
        }

        nmod_mpoly_set_str_pretty(B,
                "1 + x1*x2^2 + x2^9999999999999999999999999*x3^9", NULL, ctxB);

        nmod_mpoly_set_str_pretty(C + 0, "x1 + x2", NULL, ctxAC);
        nmod_mpoly_set_str_pretty(C + 1, "x1 - x2", NULL, ctxAC);
        nmod_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (nmod_mpoly_compose_nmod_mpoly(A, B, Cp, ctxB, ctxAC) ||
            nmod_mpoly_compose_nmod_mpoly_horner(A1, B, Cp, ctxB, ctxAC) ||
            nmod_mpoly_compose_nmod_mpoly_geobucket(A2, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 1\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_set_str_pretty(C + 0, "x1", NULL, ctxAC);
        nmod_mpoly_set_str_pretty(C + 1, "2*x2", NULL, ctxAC);
        nmod_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (!nmod_mpoly_compose_nmod_mpoly(A, B, Cp, ctxB, ctxAC) ||
            !nmod_mpoly_compose_nmod_mpoly_horner(A1, B, Cp, ctxB, ctxAC) ||
            !nmod_mpoly_compose_nmod_mpoly_geobucket(A2, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 2\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_set_str_pretty(C + 0, "2*x1", NULL, ctxAC);
        nmod_mpoly_set_str_pretty(C + 1, "x2", NULL, ctxAC);
        nmod_mpoly_set_str_pretty(C + 2, "1", NULL, ctxAC);
        if (!nmod_mpoly_compose_nmod_mpoly(A, B, Cp, ctxB, ctxAC) ||
            !nmod_mpoly_compose_nmod_mpoly_horner(A1, B, Cp, ctxB, ctxAC) ||
            !nmod_mpoly_compose_nmod_mpoly_geobucket(A2, B, Cp, ctxB, ctxAC))
        {
            printf("FAIL\n");
            flint_printf("Check example 3\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(B, ctxB);
        nmod_mpoly_clear(A, ctxAC);
        nmod_mpoly_clear(A1, ctxAC);
        nmod_mpoly_clear(A2, ctxAC);
        for (i = 0; i < 3; i++)
            nmod_mpoly_clear(C + i, ctxAC);

        nmod_mpoly_ctx_clear(ctxB);
        nmod_mpoly_ctx_clear(ctxAC);
    }

    /* check composition with generators */
    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctxB, ctxAC;
        nmod_mpoly_t B, A, A1, A2, A3;
        nmod_mpoly_struct ** C;
        slong * c;
        slong nvarsB, nvarsAC;
        slong len;
        flint_bitcnt_t exp_bits;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);

        nmod_mpoly_ctx_init_rand(ctxB, state, 20, modulus);
        nmod_mpoly_ctx_init_rand(ctxAC, state, 20, modulus);

        nvarsB = ctxB->minfo->nvars;
        nvarsAC = ctxAC->minfo->nvars;

        nmod_mpoly_init(B, ctxB);
        nmod_mpoly_init(A, ctxAC);
        nmod_mpoly_init(A1, ctxAC);
        nmod_mpoly_init(A2, ctxAC);
        nmod_mpoly_init(A3, ctxAC);

        c = (slong *) flint_malloc(nvarsB*sizeof(slong));
        C = (nmod_mpoly_struct **) flint_malloc(nvarsB * sizeof(nmod_mpoly_struct *));
        for (v = 0; v < nvarsB; v++)
        {
            C[v] = (nmod_mpoly_struct *) flint_malloc(sizeof(nmod_mpoly_struct));
            nmod_mpoly_init(C[v], ctxAC);
        }

        for (j = 0; j < 4; j++)
        {
            len = n_randint(state, 200);
            exp_bits = n_randint(state, 100) + 1;

            for (v = 0; v < nvarsB; v++)
            {
                c[v] = n_randint(state, nvarsAC + 2) - 2;
                if (c[v] >= 0)
                    nmod_mpoly_gen(C[v], c[v], ctxAC);
                else
                    nmod_mpoly_zero(C[v], ctxAC);
            }

            nmod_mpoly_randtest_bits(B, state, len, exp_bits, ctxB);

            nmod_mpoly_compose_nmod_mpoly_gen(A, B, c, ctxB, ctxAC);

            if (!nmod_mpoly_compose_nmod_mpoly(A1, B, C, ctxB, ctxAC) ||
                !nmod_mpoly_compose_nmod_mpoly_horner(A2, B, C, ctxB, ctxAC) ||
                !nmod_mpoly_compose_nmod_mpoly_geobucket(A3, B, C, ctxB, ctxAC))
            {
                printf("FAIL\n");
                flint_printf("Check composition success with generators\n"
                                                     "i: %wd, j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            if (!nmod_mpoly_equal(A, A1, ctxAC) ||
                !nmod_mpoly_equal(A, A2, ctxAC) ||
                !nmod_mpoly_equal(A, A3, ctxAC))
            {
                printf("FAIL\n");
                flint_printf("Check composition with generators\n"
                                                     "i: %wd, j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_assert_canonical(A, ctxAC);
            nmod_mpoly_assert_canonical(A1, ctxAC);
            nmod_mpoly_assert_canonical(A2, ctxAC);
            nmod_mpoly_assert_canonical(A3, ctxAC);
        }

        for (v = 0; v < nvarsB; v++)
        {
            nmod_mpoly_clear(C[v], ctxAC);
            flint_free(C[v]);
        }
        flint_free(C);
        flint_free(c);

        nmod_mpoly_clear(B, ctxB);
        nmod_mpoly_clear(A, ctxAC);
        nmod_mpoly_clear(A1, ctxAC);
        nmod_mpoly_clear(A2, ctxAC);
        nmod_mpoly_clear(A3, ctxAC);

        nmod_mpoly_ctx_clear(ctxB);
        nmod_mpoly_ctx_clear(ctxAC);
    }

    /* Check composition with identity */
    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        slong nvars, len, exp_bits;
        nmod_mpoly_struct ** vals1;
        nmod_mpoly_t f, g, g1, g2;
        nmod_mpoly_ctx_t ctx;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);
        nvars = ctx->minfo->nvars;

        vals1 = (nmod_mpoly_struct **) flint_malloc(nvars*sizeof(nmod_mpoly_struct *));
        for (v = 0; v < nvars; v++)
        {
            vals1[v] = (nmod_mpoly_struct *) flint_malloc(sizeof(nmod_mpoly_struct));
            nmod_mpoly_init(vals1[v], ctx);
            nmod_mpoly_gen(vals1[v], v, ctx);
        }

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(g1, ctx);
        nmod_mpoly_init(g2, ctx);

        len = n_randint(state, 200);
        exp_bits = n_randint(state, 300) + 1;
        nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        if (!nmod_mpoly_compose_nmod_mpoly(g, f, vals1, ctx, ctx) ||
            !nmod_mpoly_compose_nmod_mpoly_horner(g1, f, vals1, ctx, ctx) ||
            !nmod_mpoly_compose_nmod_mpoly_geobucket(g2, f, vals1, ctx, ctx) ||
            !nmod_mpoly_equal(g, g1, ctx) ||
            !nmod_mpoly_equal(g, g2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition success\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_assert_canonical(g, ctx);
        nmod_mpoly_assert_canonical(g1, ctx);
        nmod_mpoly_assert_canonical(g2, ctx);

        if (!nmod_mpoly_equal(f, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition with identity\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(g1, ctx);
        nmod_mpoly_clear(g2, ctx);

        for (v = 0; v < nvars; v++)
        {
            nmod_mpoly_clear(vals1[v], ctx);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check composition and evalall commute */
    for (i = 0; i < 40*flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx1, ctx2;
        nmod_mpoly_t f, g, g1, g2;
        nmod_mpoly_struct ** vals1;
        mp_limb_t fe, ge;
        mp_limb_t * vals2, * vals3;
        slong nvars1, nvars2;
        slong len1, len2;
        slong exp_bound1;
        flint_bitcnt_t exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx1, state, 4, modulus);
        nmod_mpoly_ctx_init_rand(ctx2, state, 8, modulus);
        nvars1 = ctx1->minfo->nvars;
        nvars2 = ctx2->minfo->nvars;

        nmod_mpoly_init(f, ctx1);
        nmod_mpoly_init(g, ctx2);
        nmod_mpoly_init(g1, ctx2);
        nmod_mpoly_init(g2, ctx2);

        len1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars1) + 1);
        len2 = n_randint(state, 10/FLINT_MAX(WORD(1), nvars2) + 2);
        exp_bound1 = n_randint(state, 12/FLINT_MAX(WORD(1), nvars1) + 1) + 2;
        exp_bits2 = n_randint(state, 200) + 1;

        nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx1);

        vals1 = (nmod_mpoly_struct **) flint_malloc(nvars1
                                                * sizeof(nmod_mpoly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (nmod_mpoly_struct *) flint_malloc(
                                                    sizeof(nmod_mpoly_struct));
            nmod_mpoly_init(vals1[v], ctx2);
            nmod_mpoly_randtest_bound(vals1[v], state, len2, exp_bits2, ctx2);
        }

        vals2 = (mp_limb_t *) flint_malloc(nvars2*sizeof(mp_limb_t));
        for (v = 0; v < nvars2; v++)
        {
            vals2[v] = n_randlimb(state);
        }

        vals3 = (mp_limb_t *) flint_malloc(nvars1*sizeof(mp_limb_t));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = nmod_mpoly_evaluate_all_ui(vals1[v], vals2, ctx2);
        }

        if (!nmod_mpoly_compose_nmod_mpoly(g, f, vals1, ctx1, ctx2) ||
            !nmod_mpoly_compose_nmod_mpoly_horner(g1, f, vals1, ctx1, ctx2) ||
            !nmod_mpoly_compose_nmod_mpoly_geobucket(g2, f, vals1, ctx1, ctx2) ||
            !nmod_mpoly_equal(g, g1, ctx2) ||
            !nmod_mpoly_equal(g, g2, ctx2))
        {
            printf("FAIL\n");
            flint_printf("Check composition success\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_assert_canonical(g, ctx2);
        nmod_mpoly_assert_canonical(g1, ctx2);
        nmod_mpoly_assert_canonical(g2, ctx2);

        fe = nmod_mpoly_evaluate_all_ui(f, vals3, ctx1);
        ge = nmod_mpoly_evaluate_all_ui(g, vals2, ctx2);

        if (fe != ge)
        {
            printf("FAIL\n");
            flint_printf("Check composition and evalall commute\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        for (v = 0; v < nvars1; v++)
        {
            nmod_mpoly_clear(vals1[v], ctx2);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        flint_free(vals2);

        flint_free(vals3);

        nmod_mpoly_clear(f, ctx1);
        nmod_mpoly_clear(g, ctx2);
        nmod_mpoly_clear(g1, ctx2);
        nmod_mpoly_clear(g2, ctx2);

        nmod_mpoly_ctx_clear(ctx1);
        nmod_mpoly_ctx_clear(ctx2);
    }

    /* Check composition with constants matches evalall */
    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx1, ctx2;
        nmod_mpoly_t f, g, g1, g2;
        nmod_mpoly_struct ** vals1;
        mp_limb_t * vals2;
        slong nvars1;
        slong len1;
        flint_bitcnt_t exp_bits1;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx1, state, 20, modulus);
        nmod_mpoly_ctx_init_rand(ctx2, state, 10, modulus);
        nvars1 = ctx1->minfo->nvars;

        nmod_mpoly_init(f, ctx1);
        nmod_mpoly_init(g, ctx2);
        nmod_mpoly_init(g1, ctx2);
        nmod_mpoly_init(g2, ctx2);

        len1 = n_randint(state, 100);
        exp_bits1 = n_randint(state, 200) + 1;

        nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx1);

        vals1 = (nmod_mpoly_struct **) flint_malloc(nvars1
                                                * sizeof(nmod_mpoly_struct *));
        vals2 = (mp_limb_t *) flint_malloc(nvars1*sizeof(mp_limb_t));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (nmod_mpoly_struct *) flint_malloc(
                                                    sizeof(nmod_mpoly_struct));
            nmod_mpoly_init(vals1[v], ctx2);
            vals2[v] = n_randlimb(state);
            nmod_mpoly_set_ui(vals1[v], vals2[v], ctx2);
        }

        if (!nmod_mpoly_compose_nmod_mpoly(g, f, vals1, ctx1, ctx2) ||
            !nmod_mpoly_compose_nmod_mpoly_geobucket(g1, f, vals1, ctx1, ctx2) ||
            !nmod_mpoly_compose_nmod_mpoly_horner(g2, f, vals1, ctx1, ctx2))
        {
            printf("FAIL\n");
            flint_printf("Check composition success\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_assert_canonical(g, ctx2);
        nmod_mpoly_assert_canonical(g1, ctx2);
        nmod_mpoly_assert_canonical(g2, ctx2);

        if (!nmod_mpoly_is_ui(g, ctx2) ||
            nmod_mpoly_get_ui(g, ctx2) != nmod_mpoly_evaluate_all_ui(f, vals2, ctx1))
        {
            printf("FAIL\n");
            flint_printf("Check composition with constants matches evalall\ni: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        for (v = 0; v < nvars1; v++)
        {
            nmod_mpoly_clear(vals1[v], ctx2);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        flint_free(vals2);

        nmod_mpoly_clear(f, ctx1);
        nmod_mpoly_clear(g, ctx2);
        nmod_mpoly_clear(g1, ctx2);
        nmod_mpoly_clear(g2, ctx2);

        nmod_mpoly_ctx_clear(ctx1);
        nmod_mpoly_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
