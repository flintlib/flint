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
    slong i, j, v;
    FLINT_TEST_INIT(state);

    flint_printf("evaluate....");
    fflush(stdout);


    /* Check repeated evalone matches evalall */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        ordering_t ord;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        fmpz_t fe;
        fmpz ** vals;
        slong * perm;
        slong nvars, len1, exp_bound1, coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_init(fe);

        perm = (slong *) flint_malloc(nvars*sizeof(slong));

        len1 = n_randint(state, 50);
        exp_bound1 = n_randint(state, 10) + 1;
        coeff_bits = n_randint(state, 100) + 1;

        vals = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            fmpz_init(vals[v]);
            fmpz_randbits(vals[v], state, 10);
            perm[v] = v;
        }

        for (j = 0; j < 2*nvars; j++)
        {
            slong a, b, c;
            a = n_randint(state, nvars);
            b = n_randint(state, nvars);
            c = perm[a];
            perm[a] = perm[b];
            perm[b] = c;
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_evaluate_all_tree_fmpz(fe, f, vals, ctx);

            for (v = 0; v < nvars; v++)
            {
                fmpz_mpoly_evaluate_one_fmpz(f, f, perm[v], vals[perm[v]], ctx);
                fmpz_mpoly_assert_canonical(f, ctx);
            }
            if (!fmpz_mpoly_equal_fmpz(f, fe, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check repeated evalone matches evalall\ni: %wd  j: %wd\n", i, j);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpz_mpoly_clear(f, ctx);

        fmpz_clear(fe);

        flint_free(perm);
    }


    /* Check addition commutes with evalall */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, fg;
        fmpz_t fe, ge, fge, t;
        ordering_t ord;
        fmpz ** vals;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(fg, ctx);

        fmpz_init(fe);
        fmpz_init(ge);
        fmpz_init(fge);
        fmpz_init(t);

        len1 = n_randint(state, 500);
        len2 = n_randint(state, 500);

        exp_bound1 = n_randint(state, 5000/nvars/nvars) + 1;
        exp_bound2 = n_randint(state, 5000/nvars/nvars) + 1;

        coeff_bits = n_randint(state, 100);


        vals = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            fmpz_init(vals[v]);
            fmpz_randbits(vals[v], state, 10);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_add(fg, f, g, ctx);

            fmpz_mpoly_evaluate_all_tree_fmpz(fe, f, vals, ctx);
            fmpz_mpoly_evaluate_all_tree_fmpz(ge, g, vals, ctx);
            fmpz_mpoly_evaluate_all_tree_fmpz(fge, fg, vals, ctx);

            fmpz_add(t, fe, ge);
            if (!fmpz_equal(t, fge))
            {
                printf("FAIL\n");
                flint_printf("Check addition commutes with evalall\ni: %wd  j: %wd\n", i, j);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(fg, ctx);

        fmpz_clear(fe);
        fmpz_clear(ge);
        fmpz_clear(fge);
        fmpz_clear(t);

    }

    /* Check multiplication commutes with evalall */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, fg;
        fmpz_t fe, ge, fge, t;
        ordering_t ord;
        fmpz ** vals;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(fg, ctx);

        fmpz_init(fe);
        fmpz_init(ge);
        fmpz_init(fge);
        fmpz_init(t);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bound1 = n_randint(state, 1000/nvars/nvars) + 1;
        exp_bound2 = n_randint(state, 1000/nvars/nvars) + 1;

        coeff_bits = n_randint(state, 100);

        vals = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            fmpz_init(vals[v]);
            fmpz_randbits(vals[v], state, 10);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_mul_johnson(fg, f, g, ctx);

            fmpz_mpoly_evaluate_all_tree_fmpz(fe, f, vals, ctx);
            fmpz_mpoly_evaluate_all_tree_fmpz(ge, g, vals, ctx);
            fmpz_mpoly_evaluate_all_tree_fmpz(fge, fg, vals, ctx);

            fmpz_mul(t, fe, ge);
            if (!fmpz_equal(t, fge))
            {
                printf("FAIL\n");
                flint_printf("Check multiplication commutes with evalall\ni: %wd  j: %wd\n", i, j);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(fg, ctx);

        fmpz_clear(fe);
        fmpz_clear(ge);
        fmpz_clear(fge);
        fmpz_clear(t);

    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

