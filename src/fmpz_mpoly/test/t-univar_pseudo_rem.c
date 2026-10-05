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

TEST_FUNCTION_START(fmpz_mpoly_univar_pseudo_rem, state)
{
    slong i, j;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, r, lc, t;
        fmpz_mpoly_univar_t ax, bx, rx;
        slong var, dega, degb, delta;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        if (ctx->minfo->nvars < 1)
        {
            fmpz_mpoly_ctx_clear(ctx);
            continue;
        }
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(r, ctx);
        fmpz_mpoly_init(lc, ctx);
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_univar_init(ax, ctx);
        fmpz_mpoly_univar_init(bx, ctx);
        fmpz_mpoly_univar_init(rx, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(a, state, 1 + n_randint(state, 10), 1 + n_randint(state, 50), 1 + n_randint(state, 4), ctx);
            do {
                fmpz_mpoly_randtest_bound(b, state, 1 + n_randint(state, 10), 1 + n_randint(state, 50), 1 + n_randint(state, 4), ctx);
            } while (fmpz_mpoly_is_zero(b, ctx));

            var = n_randint(state, ctx->minfo->nvars);
            fmpz_mpoly_to_univar(ax, a, var, ctx);
            fmpz_mpoly_to_univar(bx, b, var, ctx);
            dega = fmpz_mpoly_degree_si(a, var, ctx);
            degb = fmpz_mpoly_degree_si(b, var, ctx);

            fmpz_mpoly_univar_pseudo_rem(rx, ax, bx, ctx);
            fmpz_mpoly_from_univar(r, rx, var, ctx);

            if (dega < degb)
            {
                if (!fmpz_mpoly_equal(r, a, ctx))
                {
                    flint_printf("FAIL: deg a < deg b\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                delta = dega - degb + 1;
                fmpz_mpoly_univar_get_term_coeff(lc, bx, 0, ctx);
                fmpz_mpoly_pow_ui(lc, lc, delta, ctx);
                fmpz_mpoly_mul(t, lc, a, ctx);
                fmpz_mpoly_sub(t, t, r, ctx);

                if (fmpz_mpoly_degree_si(r, var, ctx) >= degb ||
                    !fmpz_mpoly_divides(t, t, b, ctx))
                {
                    flint_printf("FAIL: lc^delta a - r not divisible by b\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mpoly_univar_clear(ax, ctx);
        fmpz_mpoly_univar_clear(bx, ctx);
        fmpz_mpoly_univar_clear(rx, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_mpoly_clear(lc, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
