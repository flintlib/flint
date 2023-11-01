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

TEST_FUNCTION_START(nmod_mpoly_compose_nmod_poly, state)
{
    slong i, v;

    {
        nmod_poly_t A;
        nmod_mpoly_t B;
        nmod_poly_struct * Cp[3];
        nmod_poly_struct C[3];
        nmod_mpoly_ctx_t ctxB;

        nmod_mpoly_ctx_init(ctxB, 3, ORD_LEX, 13);

        nmod_mpoly_init(B, ctxB);
        nmod_poly_init_mod(A, ctxB->mod);
        for (i = 0; i < 3; i++)
        {
            Cp[i] = C + i;
            nmod_poly_init_mod(C + i, ctxB->mod);
        }

        nmod_mpoly_set_str_pretty(B,
                "1 + x1*x2^2 + x2^9999999999999999999999999*x3^9", NULL, ctxB);

        nmod_poly_zero(C + 0);
        nmod_poly_zero(C + 1);
        nmod_poly_zero(C + 2);
        nmod_poly_set_coeff_ui(C + 0, 1, 1);
        nmod_poly_set_coeff_ui(C + 1, 2, 2);
        nmod_poly_set_coeff_ui(C + 2, 3, 3);
        if (nmod_mpoly_compose_nmod_poly(A, B, Cp, ctxB))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 1\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_zero(C + 0);
        nmod_poly_zero(C + 1);
        nmod_poly_zero(C + 2);
        nmod_poly_set_coeff_ui(C + 0, 0, 1);
        nmod_poly_set_coeff_ui(C + 1, 0, 2);
        nmod_poly_set_coeff_ui(C + 2, 0, 3);
        if (!nmod_mpoly_compose_nmod_poly(A, B, Cp, ctxB))
        {
            printf("FAIL\n");
            flint_printf("Check example 2\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (!nmod_poly_is_zero(A))
        {
            printf("FAIL\n");
            flint_printf("Check example 2 equality\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(B, ctxB);
        nmod_poly_clear(A);
        for (i = 0; i < 3; i++)
            nmod_poly_clear(C + i);

        nmod_mpoly_ctx_clear(ctxB);
    }

    /* Check composition and evalall commute */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx1;
        nmod_mpoly_t f;
        nmod_poly_t g;
        nmod_poly_struct ** vals1;
        mp_limb_t fe, ge;
        mp_limb_t vals2, * vals3;
        slong nvars1;
        slong len1, len2;
        slong exp_bound1;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx1, state, 3, modulus);
        nvars1 = ctx1->minfo->nvars;

        nmod_mpoly_init(f, ctx1);
        nmod_poly_init(g, modulus);

        len1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars1) + 1);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 200/FLINT_MAX(WORD(1), nvars1) + 2) + 1;

        nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx1);

        vals1 = (nmod_poly_struct **) flint_malloc(nvars1
                                                * sizeof(nmod_poly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (nmod_poly_struct *) flint_malloc(
                                                    sizeof(nmod_poly_struct));
            nmod_poly_init(vals1[v], modulus);
            nmod_poly_randtest(vals1[v], state, len2);
        }

        vals2 = n_randint(state, modulus);

        vals3 = (mp_limb_t *) flint_malloc(nvars1*sizeof(mp_limb_t));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = nmod_poly_evaluate_nmod(vals1[v], vals2);
        }

        if (nmod_mpoly_total_degree_si(f, ctx1) < 100)
        {
            if (!nmod_mpoly_compose_nmod_poly(g, f, vals1, ctx1))
            {
                printf("FAIL\n");
                flint_printf("Check composition success\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }

            fe = nmod_mpoly_evaluate_all_ui(f, vals3, ctx1);
            ge = nmod_poly_evaluate_nmod(g, vals2);

            if (fe != ge)
            {
                printf("FAIL\n");
                flint_printf("Check composition and evalall commute\ni: %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars1; v++)
        {
            nmod_poly_clear(vals1[v]);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        flint_free(vals3);

        nmod_mpoly_clear(f, ctx1);
        nmod_poly_clear(g);

        nmod_mpoly_ctx_clear(ctx1);
    }

    TEST_FUNCTION_END(state);
}
