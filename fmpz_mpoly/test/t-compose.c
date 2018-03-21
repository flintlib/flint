/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"

int
main(void)
{
    slong i, v;
    FLINT_TEST_INIT(state);

    flint_printf("compose....");
    fflush(stdout);


    /* Check composition with identity */
    for (i = 0; i < 20*flint_test_multiplier(); i++)
    {
        slong nvars, len1, exp_bound1, coeff_bits;
        fmpz_mpoly_struct ** vals1;
        fmpz_mpoly_t f, g;
        fmpz_mpoly_ctx_t ctx;

        nvars = n_randint(state, 10) + 1;
        fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);

        vals1 = (fmpz_mpoly_struct **) flint_malloc(nvars*sizeof(fmpz_mpoly_struct *));
        for (v = 0; v < nvars; v++)
        {
            vals1[v] = (fmpz_mpoly_struct *) flint_malloc(sizeof(fmpz_mpoly_struct)); 
            fmpz_mpoly_init(vals1[v], ctx);
            fmpz_mpoly_gen(vals1[v], v, ctx);            
        }

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);

        len1 = n_randint(state, 200);
        exp_bound1 = n_randint(state, 40) + 1;
        coeff_bits = n_randint(state, 100) + 1;
        fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
        fmpz_mpoly_compose(g, f, vals1, ctx, ctx);

        if (!fmpz_mpoly_equal(f, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check composition with identity\ni: %wd\n", i);
            flint_abort();
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(f, ctx);

        for (v = 0; v < nvars; v++)
        {
            fmpz_mpoly_clear(vals1[v], ctx);
            flint_free(vals1[v]);            
        }
        flint_free(vals1);
    }

    /* Check composition and evalall commute */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        ordering_t ord1, ord2;
        fmpz_mpoly_ctx_t ctx1, ctx2;
        fmpz_mpoly_t f, g;
        fmpz_mpoly_struct ** vals1;
        fmpz_t fe, ge;
        fmpz ** vals2, ** vals3;
        slong nvars1, nvars2;
        slong len1, len2;
        slong exp_bound1, exp_bound2;
        slong coeff_bits;

        ord1 = mpoly_ordering_randtest(state);
        ord2 = mpoly_ordering_randtest(state);
        nvars1 = n_randint(state, 10) + 1;
        nvars2 = n_randint(state, 10) + 1;
        fmpz_mpoly_ctx_init(ctx1, nvars1, ord1);
        fmpz_mpoly_ctx_init(ctx2, nvars2, ord2);

        fmpz_mpoly_init(f, ctx1);
        fmpz_mpoly_init(g, ctx2);
        fmpz_init(fe);
        fmpz_init(ge);

        len1 = n_randint(state, 80/nvars1 + 1);
        len2 = n_randint(state, 30/nvars2 + 1);
        exp_bound1 = n_randint(state, 15/nvars1 + 2) + 1;
        exp_bound2 = n_randint(state, 15/nvars2 + 2) + 1;
        coeff_bits = n_randint(state, 10);

        vals1 = (fmpz_mpoly_struct **) flint_malloc(nvars1
                                                * sizeof(fmpz_mpoly_struct *));
        for (v = 0; v < nvars1; v++)
        {
            vals1[v] = (fmpz_mpoly_struct *) flint_malloc(
                                                    sizeof(fmpz_mpoly_struct)); 
            fmpz_mpoly_init(vals1[v], ctx2);
            fmpz_mpoly_randtest_bound(vals1[v], state, len2,
                                                 coeff_bits, exp_bound2, ctx2);
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
            fmpz_mpoly_evaluate_all_tree_fmpz(vals3[v], vals1[v], vals2, ctx2);
        }

        fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx1);

        fmpz_mpoly_compose(g, f, vals1, ctx1, ctx2);
        fmpz_mpoly_assert_canonical(g, ctx2);

        fmpz_mpoly_evaluate_all_tree_fmpz(fe, f, vals3, ctx1);
        fmpz_mpoly_evaluate_all_tree_fmpz(ge, g, vals2, ctx2);

        if (!fmpz_equal(fe, ge))
        {
            printf("FAIL\n");
            flint_printf("Check composition and evalall commute\ni: %wd\n", i);
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

        fmpz_clear(fe);
        fmpz_clear(ge);

    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

