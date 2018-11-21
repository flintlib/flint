/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"

int
main(void)
{
    slong i, v;
    FLINT_TEST_INIT(state);

    flint_printf("compose_nmod_mpoly....");
    fflush(stdout);

    /* Check composition with identity */
    for (i = 0; i < 20*flint_test_multiplier(); i++)
    {
        slong nvars, len, exp_bits;
        nmod_mpoly_struct ** vals1;
        nmod_mpoly_t f, g;
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

        len = n_randint(state, 200);
        exp_bits = n_randint(state, 300) + 1;
        nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
        nmod_mpoly_compose_nmod_mpoly(g, f, vals1, ctx, ctx);
        nmod_mpoly_assert_canonical(g, ctx);

        if (!nmod_mpoly_equal(f, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition with identity\ni: %wd\n", i);
            flint_abort();
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(f, ctx);

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
        nmod_mpoly_t f, g;
        nmod_mpoly_struct ** vals1;
        mp_limb_t fe, ge;
        mp_limb_t * vals2, * vals3;
        slong nvars1, nvars2;
        slong len1, len2;
        slong exp_bound1;
        mp_bitcnt_t exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx1, state, 4, modulus);
        nmod_mpoly_ctx_init_rand(ctx2, state, 8, modulus);
        nvars1 = ctx1->minfo->nvars;
        nvars2 = ctx2->minfo->nvars;

        nmod_mpoly_init(f, ctx1);
        nmod_mpoly_init(g, ctx2);

        len1 = n_randint(state, 50/nvars1 + 1);
        len2 = n_randint(state, 20/nvars2 + 1);
        exp_bound1 = n_randint(state, 15/nvars1 + 2) + 1;
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

        if (nmod_mpoly_total_degree_si(f, ctx1)
                                         <= 4000/(1+len2*len2*nvars2*nvars2))
        {
            nmod_mpoly_compose_nmod_mpoly(g, f, vals1, ctx1, ctx2);
            nmod_mpoly_assert_canonical(g, ctx2);

            fe = nmod_mpoly_evaluate_all_ui(f, vals3, ctx1);
            ge = nmod_mpoly_evaluate_all_ui(g, vals2, ctx2);

            if (fe != ge)
            {
                printf("FAIL\n");
                flint_printf("Check composition and evalall commute\ni: %wd\n", i);
                flint_abort();
            }
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

        nmod_mpoly_ctx_clear(ctx1);
        nmod_mpoly_ctx_clear(ctx2);
    }

    /* Check composition with constants matches evalall */
    for (i = 0; i < 40*flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx1, ctx2;
        nmod_mpoly_t f, g;
        nmod_mpoly_struct ** vals1;
        mp_limb_t * vals2;
        slong nvars1;
        slong len1;
        mp_bitcnt_t exp_bits1;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx1, state, 20, modulus);
        nmod_mpoly_ctx_init_rand(ctx2, state, 10, modulus);
        nvars1 = ctx1->minfo->nvars;

        nmod_mpoly_init(f, ctx1);
        nmod_mpoly_init(g, ctx2);

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

        nmod_mpoly_compose_nmod_mpoly(g, f, vals1, ctx1, ctx2);
        nmod_mpoly_assert_canonical(g, ctx2);

        if (!nmod_mpoly_is_ui(g, ctx2)
            || nmod_mpoly_get_ui(g, ctx2) 
                    != nmod_mpoly_evaluate_all_ui(f, vals2, ctx1))
        {
            printf("FAIL\n");
            flint_printf("Check composition with constants matches evalall\ni: %wd\n", i);
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

        nmod_mpoly_ctx_clear(ctx1);
        nmod_mpoly_ctx_clear(ctx2);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}


