/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "calcium.h"
#include "utils_flint.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("fmpz_mpoly_set_gen_fmpz_poly...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t F, G, x, t;
        fmpz_poly_t f;
        slong nvars, i, j;

        fmpz_poly_init(f);

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_init(F, ctx);
        fmpz_mpoly_init(G, ctx);
        fmpz_mpoly_init(x, ctx);
        fmpz_mpoly_init(t, ctx);

        fmpz_mpoly_randtest_bound(F, state, 10, 100, 100, ctx);
        fmpz_poly_randtest(f, state, 10, 100);

        i = n_randint(state, nvars);
        fmpz_mpoly_set_gen_fmpz_poly(F, i, f, ctx);

        fmpz_mpoly_gen(x, i, ctx);
        for (j = 0; j < f->length; j++)
        {
            fmpz_mpoly_pow_ui(t, x, j, ctx);
            fmpz_mpoly_scalar_mul_fmpz(t, t, f->coeffs + j, ctx);
            fmpz_mpoly_add(G, G, t, ctx);
        }

        if (!fmpz_mpoly_equal(F, G, ctx))
        {
            flint_printf("FAIL\n");
            fmpz_mpoly_print_pretty(F, NULL, ctx); flint_printf("\n\n");
            fmpz_mpoly_print_pretty(G, NULL, ctx); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_mpoly_clear(F, ctx);
        fmpz_mpoly_clear(G, ctx);
        fmpz_mpoly_clear(x, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
