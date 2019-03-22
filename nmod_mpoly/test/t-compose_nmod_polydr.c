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

    flint_printf("compose_nmod_poly....");
    fflush(stdout);

    /* Check composition and evalall commute */
    for (i = 0; i < 50*flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx1;
        nmod_mpoly_t f;
        nmod_polydr_t g;
        nmod_polydr_struct ** vals1;
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
        nmod_polydr_init(g, ctx1->ffinfo);

        len1 = n_randint(state, 50/nvars1 + 1);
        len2 = n_randint(state, 100);
        exp_bound1 = n_randint(state, 200/nvars1 + 2) + 1;

        nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx1);

        vals1 = (nmod_polydr_struct **) flint_malloc(nvars1
                                                * sizeof(nmod_polydr_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (nmod_polydr_struct *) flint_malloc(
                                                    sizeof(nmod_polydr_struct)); 
            nmod_polydr_init(vals1[v], ctx1->ffinfo);
            nmod_polydr_randtest(vals1[v], state, len2, ctx1->ffinfo);
        }

        vals2 = n_randint(state, modulus);

        vals3 = (mp_limb_t *) flint_malloc(nvars1*sizeof(mp_limb_t));
        for (v = 0; v < nvars1; v++)
        {
            vals3[v] = nmod_polydr_evaluate_nmod(vals1[v], vals2, ctx1->ffinfo);
        }

        if (nmod_mpoly_total_degree_si(f, ctx1) < 100)
        {
            nmod_mpoly_compose_nmod_polydr(g, f, vals1, ctx1);

            fe = nmod_mpoly_evaluate_all_ui(f, vals3, ctx1);
            ge = nmod_polydr_evaluate_nmod(g, vals2, ctx1->ffinfo);

            if (fe != ge)
            {
                printf("FAIL\n");
                flint_printf("Check composition and evalall commute\ni: %wd\n", i);
                flint_abort();
            }
        }

        for (v = 0; v < nvars1; v++)
        {
            nmod_polydr_clear(vals1[v], ctx1->ffinfo);
            flint_free(vals1[v]);
        }
        flint_free(vals1);

        flint_free(vals3);

        nmod_mpoly_clear(f, ctx1);
        nmod_polydr_clear(g, ctx1->ffinfo);

        nmod_mpoly_ctx_clear(ctx1);
    }
    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
