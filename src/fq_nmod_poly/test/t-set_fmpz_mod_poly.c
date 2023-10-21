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
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fq_nmod.h"
#include "fq_nmod_poly.h"

TEST_FUNCTION_START(fq_nmod_poly_set_fmpz_mod_poly, state)
{
    int i, result;

    /* Check litfed polynomials by evaluating at random points */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong len;
        fq_nmod_ctx_t ctx;
        fq_nmod_t r, s;
        fq_nmod_poly_t a;
        fmpz_mod_poly_t b;
        fmpz_t p;
        fmpz_mod_ctx_t ctxp;

        len = n_randint(state, 15) + 1;
        fq_nmod_ctx_randtest(ctx, state);
        fq_nmod_init(r, ctx);
        fq_nmod_init(s, ctx);
        fq_nmod_poly_init(a, ctx);
        fmpz_mod_ctx_init(ctxp, fq_nmod_ctx_prime(ctx));
        fmpz_mod_poly_init(b, ctxp);
        fmpz_init(p);

        fmpz_mod_poly_randtest(b, state, len, ctxp);
        fmpz_randtest(p, state, 10);

        fq_nmod_poly_set_fmpz_mod_poly(a, b, ctx);
        fq_nmod_set_fmpz(r, p, ctx);
        fq_nmod_poly_evaluate_fq_nmod(r, a, r, ctx);
        fmpz_mod_poly_evaluate_fmpz(p, b, p, ctxp);
        fq_nmod_set_fmpz(s, p, ctx);

        result = fq_nmod_equal(r, s, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"); fq_nmod_ctx_print(ctx);
            flint_printf("\nb = "); fmpz_mod_poly_print_pretty(b, "X", ctxp);
            flint_printf("\np = "); fmpz_print(p); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(r, ctx);
        fq_nmod_clear(s, ctx);
        fmpz_mod_poly_clear(b, ctxp);
        fmpz_mod_ctx_clear(ctxp);
        fmpz_clear(p);
        fq_nmod_poly_clear(a, ctx);
        fq_nmod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
