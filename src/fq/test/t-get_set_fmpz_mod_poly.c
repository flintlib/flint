/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

TEST_FUNCTION_START(fq_get_set_fmpz_mod_poly, state)
{
    slong i, j;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_t x, y;
        fmpz_mod_poly_t z, t;

        fq_ctx_randtest(ctx, state);
        fq_init(x, ctx);
        fq_init(y, ctx);
        fmpz_mod_poly_init(z, ctx->ctxp);
        fmpz_mod_poly_init(t, ctx->ctxp);

        for (j = 0; j < 100; j++)
        {
            fq_rand(x, state, ctx);
            fq_rand(y, state, ctx);
            fq_get_fmpz_mod_poly(z, x, ctx);
            fmpz_mod_poly_randtest(t, state, 20, ctx->ctxp);
            fmpz_mod_poly_mul(t, t, ctx->modulus, ctx->ctxp);
            fmpz_mod_poly_add(z, z, t, ctx->ctxp);
            fq_set_fmpz_mod_poly(y, z, ctx);

            if (!fq_equal(y, x, ctx))
            {
                flint_printf("FAIL:\n");
                flint_printf("check get/set match i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_poly_clear(t, ctx->ctxp);
        fmpz_mod_poly_clear(z, ctx->ctxp);
        fq_clear(y, ctx);
        fq_clear(x, ctx);
        fq_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
