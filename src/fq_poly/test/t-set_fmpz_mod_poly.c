/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_poly.h"

TEST_FUNCTION_START(fq_poly_set_fmpz_mod_poly, state)
{
    int i, result;

    /* Check litfed polynomials by evaluating at random points */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong len;
        fq_ctx_t ctx;
        fq_t r, s;
        fq_poly_t a;
        fmpz_mod_poly_t b;
        fmpz_t p;

        len = n_randint(state, 15) + 1;
        fq_ctx_randtest(ctx, state);
        fq_init(r, ctx);
        fq_init(s, ctx);
        fq_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx->ctxp);
        fmpz_init(p);

        fmpz_mod_poly_randtest(b, state, len, ctx->ctxp);
        fmpz_randtest(p, state, 10);

        fq_poly_set_fmpz_mod_poly(a, b, ctx);
        fq_set_fmpz(r, p, ctx);
        fq_poly_evaluate_fq(r, a, r, ctx);
        fmpz_mod_poly_evaluate_fmpz(p, b, p, ctx->ctxp);
        fq_set_fmpz(s, p, ctx);

        result = fq_equal(r, s, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n");
            fq_ctx_print(ctx);
            flint_printf("\nb = "); fmpz_mod_poly_print_pretty(b, "X", ctx->ctxp);
            flint_printf("\np = "); fmpz_print(p); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fq_clear(r, ctx);
        fq_clear(s, ctx);
        fmpz_mod_poly_clear(b, ctx->ctxp);
        fmpz_clear(p);
        fq_poly_clear(a, ctx);
        fq_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
