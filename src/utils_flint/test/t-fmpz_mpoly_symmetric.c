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

    flint_printf("fmpz_mpoly_symmetric...");
    fflush(stdout);

    flint_randinit(state);

    /* todo; this does noting useful right now except verifying
       that the function can be called */
    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t F;
        slong nvars, n, k;
        slong vars[10];

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        nvars = ctx->minfo->nvars;

        n = n_randint(state, nvars + 1);
        for (k = 0; k < n; k++)
            vars[k] = k;

        k = n_randint(state, 5);
        fmpz_mpoly_init(F, ctx);
        fmpz_mpoly_symmetric_gens(F, k, vars, n, ctx);

        fmpz_mpoly_clear(F, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
