/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_compose_fmpz_poly, state)
{
    slong i, v;

    {
        fmpz_poly_t A;
        fmpz_mpoly_t B;
        fmpz_poly_struct * Cp[3];
        fmpz_poly_struct C[3];
        fmpz_mpoly_ctx_t ctxB;

        fmpz_mpoly_ctx_init(ctxB, 3, ORD_LEX);

        fmpz_mpoly_init(B, ctxB);
        fmpz_poly_init(A);
        for (i = 0; i < 3; i++)
        {
            Cp[i] = C + i;
            fmpz_poly_init(C + i);
        }

        fmpz_mpoly_set_str_pretty(B,
                "1 + x1*x2^2 + x2^9999999999999999999999999*x3^9", NULL, ctxB);

        fmpz_poly_zero(C + 0);
        fmpz_poly_zero(C + 1);
        fmpz_poly_zero(C + 2);
        fmpz_poly_set_coeff_si(C + 0, 1, 1);
        fmpz_poly_set_coeff_si(C + 1, 2, 2);
        fmpz_poly_set_coeff_si(C + 2, 3, 3);
        if (fmpz_mpoly_compose_fmpz_poly(A, B, Cp, ctxB))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 1\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_zero(C + 0);
        fmpz_poly_zero(C + 1);
        fmpz_poly_zero(C + 2);
        fmpz_poly_set_coeff_si(C + 0, 0, 1);
        fmpz_poly_set_coeff_si(C + 1, 0, 1);
        fmpz_poly_set_coeff_si(C + 2, 0, 1);
        if (!fmpz_mpoly_compose_fmpz_poly(A, B, Cp, ctxB))
        {
            printf("FAIL\n");
            flint_printf("Check example 2\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_zero(C + 0);
        fmpz_poly_set_coeff_si(C + 0, 0, 3);
        if (!fmpz_poly_equal(A, C + 0))
        {
            printf("FAIL\n");
            flint_printf("Check example 2 equality\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(B, ctxB);
        fmpz_poly_clear(A);
        for (i = 0; i < 3; i++)
            fmpz_poly_clear(C + i);

        fmpz_mpoly_ctx_clear(ctxB);
    }

    /* Check composition and evalall commute */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx1;
        fmpz_mpoly_t f;
        fmpz_poly_t g;
        fmpz_poly_struct ** vals1;
        fmpz_t fe, ge;
        fmpz_t vals2;
        fmpz ** vals3;
        slong nvars1;
        slong len1, len2;
        slong exp_bound1;
        flint_bitcnt_t coeff_bits, coeff_bits2;

        fmpz_mpoly_ctx_init_rand(ctx1, state, 10);
        nvars1 = ctx1->minfo->nvars;

        fmpz_mpoly_init(f, ctx1);
        fmpz_poly_init(g);

        fmpz_init(fe);
        fmpz_init(ge);

        len1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars1) + 1);
        len2 = n_randint(state, 40);
        exp_bound1 = n_randint(state, 200/FLINT_MAX(WORD(1), nvars1) + 2) + 1;
        coeff_bits = n_randint(state, 100) + 1;
        coeff_bits2 = n_randint(state, 10) + 1;

        fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx1);

        vals1 = (fmpz_poly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpz_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpz_poly_struct *) flint_malloc(
                                                    sizeof(fmpz_poly_struct));
            fmpz_poly_init(vals1[v]);
            fmpz_poly_randtest(vals1[v], state, len2, coeff_bits2);
        }

        fmpz_init(vals2);
        fmpz_randbits(vals2, state, 10);

        vals3 = (fmpz **) flint_malloc(nvars1*sizeof(fmpz *));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals3[v]);
            fmpz_poly_evaluate_fmpz(vals3[v], vals1[v], vals2);
        }

        if (fmpz_mpoly_total_degree_si(f, ctx1) < 50)
        {
            if (!fmpz_mpoly_compose_fmpz_poly(g, f, vals1, ctx1) ||
                !fmpz_mpoly_evaluate_all_fmpz(fe, f, vals3, ctx1))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_evaluate_fmpz(ge, g, vals2);

            if (!fmpz_equal(fe, ge))
            {
                printf("FAIL\n");
                flint_printf("Check composition and evalall commute\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpz_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        fmpz_clear(vals2);

        for (v = 0; v < nvars1; v++)
        {
            fmpz_poly_clear(vals1[v]);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fmpz_clear(fe);
        fmpz_clear(ge);

        fmpz_mpoly_clear(f, ctx1);
        fmpz_poly_clear(g);

        fmpz_mpoly_ctx_clear(ctx1);
    }

    /* Check composition with constants matches evalall */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx1;
        fmpz_mpoly_t f;
        fmpz_poly_t g;
        fmpz_poly_struct ** vals1;
        fmpz_t fe;
        fmpz ** vals3;
        slong nvars1;
        slong len1;
        slong exp_bound1;
        flint_bitcnt_t coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx1, state, 10);
        nvars1 = ctx1->minfo->nvars;

        fmpz_mpoly_init(f, ctx1);
        fmpz_poly_init(g);

        fmpz_init(fe);

        len1 = n_randint(state, 50);
        exp_bound1 = n_randint(state, 200/FLINT_MAX(WORD(1), nvars1) + 2) + 1;
        coeff_bits = n_randint(state, 100) + 1;

        fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx1);

        vals3 = (fmpz **) flint_malloc(nvars1*sizeof(fmpz *));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals3[v]);
            fmpz_randtest(vals3[v], state, 20);
        }

        vals1 = (fmpz_poly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpz_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpz_poly_struct *) flint_malloc(
                                                    sizeof(fmpz_poly_struct));
            fmpz_poly_init(vals1[v]);
            fmpz_poly_set_fmpz(vals1[v], vals3[v]);
        }

        if (fmpz_mpoly_total_degree_si(f, ctx1) < 50)
        {
            fmpz_poly_t t;

            if (!fmpz_mpoly_compose_fmpz_poly(g, f, vals1, ctx1) ||
                !fmpz_mpoly_evaluate_all_fmpz(fe, f, vals3, ctx1))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_init(t);
            fmpz_poly_set_fmpz(t, fe);
            if (!fmpz_poly_equal(g, t))
            {
                printf("FAIL\n");
                flint_printf("Check composition with constants matches evalall\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            fmpz_poly_clear(t);
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpz_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        for (v = 0; v < nvars1; v++)
        {
            fmpz_poly_clear(vals1[v]);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fmpz_clear(fe);

        fmpz_mpoly_clear(f, ctx1);
        fmpz_poly_clear(g);

        fmpz_mpoly_ctx_clear(ctx1);
    }

    /* Check multiprecision composition with constants matches evalall */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx1;
        fmpz_mpoly_t f;
        fmpz_poly_t g;
        fmpz_poly_struct ** vals1;
        fmpz_t fe;
        fmpz ** vals3;
        slong nvars1;
        slong len1;
        flint_bitcnt_t exp_bits;
        flint_bitcnt_t coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx1, state, 10);
        nvars1 = ctx1->minfo->nvars;

        fmpz_mpoly_init(f, ctx1);
        fmpz_poly_init(g);

        fmpz_init(fe);

        len1 = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 100) + 1;

        fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx1);

        vals3 = (fmpz **) flint_malloc(nvars1*sizeof(fmpz *));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals3[v]);
            fmpz_set_si(vals3[v], n_randint(state, UWORD(3)) - WORD(1));
        }

        vals1 = (fmpz_poly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpz_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpz_poly_struct *) flint_malloc(
                                                    sizeof(fmpz_poly_struct));
            fmpz_poly_init(vals1[v]);
            fmpz_poly_set_fmpz(vals1[v], vals3[v]);
        }

        {
            fmpz_poly_t t;

            if (!fmpz_mpoly_compose_fmpz_poly(g, f, vals1, ctx1) ||
                !fmpz_mpoly_evaluate_all_fmpz(fe, f, vals3, ctx1))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_init(t);
            fmpz_poly_set_fmpz(t, fe);
            if (!fmpz_poly_equal(g, t))
            {
                printf("FAIL\n");
                flint_printf("Check multiprecision composition with constants matches evalall\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            fmpz_poly_clear(t);
        }

        for (v = 0; v < nvars1; v++)
        {
            fmpz_clear(vals3[v]);
            flint_free(vals3[v]);
        }
        flint_free(vals3);

        for (v = 0; v < nvars1; v++)
        {
            fmpz_poly_clear(vals1[v]);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        fmpz_clear(fe);

        fmpz_mpoly_clear(f, ctx1);
        fmpz_poly_clear(g);

        fmpz_mpoly_ctx_clear(ctx1);
    }

    TEST_FUNCTION_END(state);
}
