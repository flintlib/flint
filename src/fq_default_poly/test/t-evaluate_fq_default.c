/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default_poly.h"

TEST_FUNCTION_START(fq_default_poly_evaluate, state)
{
    int i, result;
    fmpz_t p;
    fq_default_ctx_t ctx;

    fmpz_init(p);

    /* Given a random fq_default_ctx, compute two random polynomials f1, f2 
       and a random element a. Ensure f1(a) * f2(b) == (f1 * f2)(a). */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_default_t a, b, c;
        fq_default_poly_t f1, f2, f3;

        /* Initalise a random fq_context */
        fq_default_ctx_init_randtest(ctx, state);

        /* Initalise fq_default elements */
        fq_default_init(a, ctx);
        fq_default_init(b, ctx);
        fq_default_init(c, ctx);
        fq_default_poly_init(f1, ctx);
        fq_default_poly_init(f2, ctx);
        fq_default_poly_init(f3, ctx);

        /* Set a to be a random element */
        fq_default_randtest(a, state, ctx);

        /* Set f1, f2 to be random polynomials and f3 their product */
        fq_default_poly_randtest(f1, state, n_randint(state, 10), ctx);
        fq_default_poly_randtest(f2, state, n_randint(state, 10), ctx);
        fq_default_poly_mul(f3, f1, f2, ctx);

        /* Evaluate the polynomials f1, f2, f3 = f1 * f2 on the element a */
        fq_default_poly_evaluate_fq_default(c, f3, a, ctx);
        fq_default_poly_evaluate_fq_default(b, f2, a, ctx);
        fq_default_poly_evaluate_fq_default(a, f1, a, ctx);

        /* Ensure that f1(a) * f2(a) = f3(a) = (f1 * f2)(a) */
        fq_default_mul(a, a, b, ctx); // a = f1(a) * f2(a)
        result = (fq_default_equal(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fq_default_ctx_print(ctx), flint_printf("\n\n");
            fq_default_poly_print(f1, ctx), flint_printf("\n\n");
            fq_default_poly_print(f2, ctx), flint_printf("\n\n");
            fq_default_poly_print(f3, ctx), flint_printf("\n\n");
            fq_default_print(a, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Clear everything */
        fq_default_poly_clear(f1, ctx);
        fq_default_poly_clear(f2, ctx);
        fq_default_poly_clear(f3, ctx);
        fq_default_clear(a, ctx);
        fq_default_clear(b, ctx);
        fq_default_clear(c, ctx);
        fq_default_ctx_clear(ctx);
    }

    fmpz_clear(p);

    TEST_FUNCTION_END(state);
}
