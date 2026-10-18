/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_univar_set_coeff_ui, state)
{
    slong i, j, k;

    /* Check that setting the coefficients of a univar in any order,
       replacing and removing some of them, gives the right univar */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, c, t;
        fmpz_mpoly_univar_t fx, gx;
        slong n, nextra, var;
        slong * perm;
        ulong * extra;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);
        if (ctx->minfo->nvars < 1)
        {
            fmpz_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_univar_init(fx, ctx);
        fmpz_mpoly_univar_init(gx, ctx);

        fmpz_mpoly_randtest_bound(f, state, n_randint(state, 30), 50, 10, ctx);
        var = n_randint(state, ctx->minfo->nvars);
        fmpz_mpoly_to_univar(fx, f, var, ctx);
        n = fmpz_mpoly_univar_length(fx, ctx);

        perm = (slong *) flint_malloc((n + 1)*sizeof(slong));
        extra = (ulong *) flint_malloc((n + 1)*sizeof(ulong));
        for (j = 0; j < n; j++)
            perm[j] = j;
        for (j = 0; j < n; j++)
        {
            k = j + n_randint(state, n - j);
            FLINT_SWAP(slong, perm[j], perm[k]);
        }

        nextra = 0;
        for (j = 0; j < n; j++)
        {
            ulong e = fmpz_mpoly_univar_get_term_exp_si(fx, perm[j], ctx);

            /* set a coefficient that is replaced below */
            if (n_randint(state, 2))
            {
                fmpz_mpoly_randtest_bound(t, state, 5, 50, 10, ctx);
                fmpz_mpoly_univar_set_coeff_ui(gx, e, t, ctx);
            }

            /* set a coefficient that is removed below */
            if (n_randint(state, 3) == 0)
            {
                extra[nextra] = 100 + n_randint(state, 100);
                fmpz_mpoly_randtest_bound(t, state, 5, 50, 10, ctx);
                fmpz_mpoly_univar_set_coeff_ui(gx, extra[nextra], t, ctx);
                nextra++;
            }

            fmpz_mpoly_univar_get_term_coeff(c, fx, perm[j], ctx);
            fmpz_mpoly_univar_set_coeff_ui(gx, e, c, ctx);
            fmpz_mpoly_univar_assert_canonical(gx, ctx);
        }

        fmpz_mpoly_zero(t, ctx);
        for (j = 0; j < nextra; j++)
            fmpz_mpoly_univar_set_coeff_ui(gx, extra[j], t, ctx);

        fmpz_mpoly_univar_assert_canonical(gx, ctx);
        fmpz_mpoly_from_univar(g, gx, var, ctx);
        if (fmpz_mpoly_univar_length(gx, ctx) != n || !fmpz_mpoly_equal(f, g, ctx))
            TEST_FUNCTION_FAIL("i = %wd, n = %wd\n", i, n);

        flint_free(perm);
        flint_free(extra);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_univar_clear(fx, ctx);
        fmpz_mpoly_univar_clear(gx, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
