/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_compose_fmpq_poly, state)
{
    slong i, v;

    /* Check composition and evalall commute */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx1;
        fmpq_mpoly_t f;
        fmpq_poly_t g;
        fmpq_poly_struct ** vals1;
        fmpq_t fe, ge;
        fmpq_t vals2;
        fmpq ** vals3;
        slong nvars1;
        slong len1, len2;
        slong exp_bound1;
        flint_bitcnt_t coeff_bits, coeff_bits2;

        fmpq_mpoly_ctx_init_rand(ctx1, state, 10);
        nvars1 = ctx1->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx1);
        fmpq_poly_init(g);

        fmpq_init(fe);
        fmpq_init(ge);

        len1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars1) + 1);
        len2 = n_randint(state, 40);
        exp_bound1 = n_randint(state, 100/FLINT_MAX(WORD(1), nvars1) + 2) + 1;
        coeff_bits = n_randint(state, 100) + 1;
        coeff_bits2 = n_randint(state, 15) + 1;

        fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx1);

        vals1 = (fmpq_poly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpq_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpq_poly_struct *) flint_malloc(
                                                    sizeof(fmpq_poly_struct));
            fmpq_poly_init(vals1[v]);
            fmpq_poly_randtest(vals1[v], state, len2, coeff_bits2);
        }

        fmpq_init(vals2);
        fmpq_randbits(vals2, state, 10);

        vals3 = (fmpq **) flint_malloc(nvars1*sizeof(fmpq *));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpq *) flint_malloc(sizeof(fmpq));
            fmpq_init(vals3[v]);
            fmpq_poly_evaluate_fmpq(vals3[v], vals1[v], vals2);
        }

        if (fmpq_mpoly_total_degree_si(f, ctx1) < 50)
        {
            if (!fmpq_mpoly_compose_fmpq_poly(g, f, vals1, ctx1) ||
                !fmpq_mpoly_evaluate_all_fmpq(fe, f, vals3, ctx1))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            fmpq_poly_evaluate_fmpq(ge, g, vals2);
            if (!fmpq_equal(fe, ge))
            {
                printf("FAIL\n");
                flint_printf("Check composition and evalall commute\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpq_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        fmpq_clear(vals2);

        for (v = 0; v < nvars1; v++)
        {
            fmpq_poly_clear(vals1[v]);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fmpq_clear(fe);
        fmpq_clear(ge);

        fmpq_mpoly_clear(f, ctx1);
        fmpq_poly_clear(g);

        fmpq_mpoly_ctx_clear(ctx1);
    }

    /* Check composition with constants matches evalall */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx1;
        fmpq_mpoly_t f;
        fmpq_poly_t g;
        fmpq_poly_struct ** vals1;
        fmpq_t fe;
        fmpq ** vals3;
        slong nvars1;
        slong len1;
        slong exp_bound1;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx1, state, 10);
        nvars1 = ctx1->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx1);
        fmpq_poly_init(g);

        fmpq_init(fe);

        len1 = n_randint(state, 50);
        exp_bound1 = n_randint(state, 200/FLINT_MAX(WORD(1), nvars1) + 2) + 1;
        coeff_bits = n_randint(state, 100) + 1;

        fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx1);

        vals3 = (fmpq **) flint_malloc(nvars1*sizeof(fmpq *));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpq *) flint_malloc(sizeof(fmpq));
            fmpq_init(vals3[v]);
            fmpq_randtest(vals3[v], state, 20);
        }

        vals1 = (fmpq_poly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpq_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpq_poly_struct *) flint_malloc(
                                                    sizeof(fmpq_poly_struct));
            fmpq_poly_init(vals1[v]);
            fmpq_poly_set_fmpq(vals1[v], vals3[v]);
        }

        if (fmpq_mpoly_total_degree_si(f, ctx1) < 50)
        {
            fmpq_poly_t t;
            if (!fmpq_mpoly_compose_fmpq_poly(g, f, vals1, ctx1) ||
                !fmpq_mpoly_evaluate_all_fmpq(fe, f, vals3, ctx1))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            fmpq_poly_init(t);
            fmpq_poly_set_fmpq(t, fe);
            if (!fmpq_poly_equal(g, t))
            {
                printf("FAIL\n");
                flint_printf("Check composition with constants matches evalall\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            fmpq_poly_clear(t);
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpq_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        for (v = 0; v < nvars1; v++)
        {
            fmpq_poly_clear(vals1[v]);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fmpq_clear(fe);

        fmpq_mpoly_clear(f, ctx1);
        fmpq_poly_clear(g);

        fmpq_mpoly_ctx_clear(ctx1);
    }

    /* Check multiprecision composition with constants matches evalall */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx1;
        fmpq_mpoly_t f;
        fmpq_poly_t g;
        fmpq_poly_struct ** vals1;
        fmpq_t fe;
        fmpq ** vals3;
        slong nvars1;
        slong len1;
        flint_bitcnt_t exp_bits;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx1, state, 10);
        nvars1 = ctx1->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx1);
        fmpq_poly_init(g);

        fmpq_init(fe);

        len1 = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 100) + 1;

        fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx1);

        vals3 = (fmpq **) flint_malloc(nvars1*sizeof(fmpq *));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpq *) flint_malloc(sizeof(fmpq));
            fmpq_init(vals3[v]);
            fmpq_set_si(vals3[v], n_randint(state, UWORD(3)) - WORD(1), UWORD(1));
        }

        vals1 = (fmpq_poly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpq_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpq_poly_struct *) flint_malloc(
                                                    sizeof(fmpq_poly_struct));
            fmpq_poly_init(vals1[v]);
            fmpq_poly_set_fmpq(vals1[v], vals3[v]);
        }

        {
            fmpq_poly_t t;
            if (!fmpq_mpoly_compose_fmpq_poly(g, f, vals1, ctx1) ||
                !fmpq_mpoly_evaluate_all_fmpq(fe, f, vals3, ctx1))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            fmpq_poly_init(t);
            fmpq_poly_set_fmpq(t, fe);
            if (!fmpq_poly_equal(g, t))
            {
                printf("FAIL\n");
                flint_printf("Check multiprecision composition with constants matches evalall\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            fmpq_poly_clear(t);
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpq_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        for (v = 0; v < nvars1; v++)
        {
            fmpq_poly_clear(vals1[v]);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fmpq_clear(fe);

        fmpq_mpoly_clear(f, ctx1);
        fmpq_poly_clear(g);

        fmpq_mpoly_ctx_clear(ctx1);
    }

    TEST_FUNCTION_END(state);
}
