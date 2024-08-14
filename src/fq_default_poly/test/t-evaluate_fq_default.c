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

// Given a context, compute a random polynomial and evaluate it twice on a random fq_element
void
test_eval_with_ctx(fq_default_ctx_t ctx,
                   flint_rand_t state)
{
    int i, result;
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_default_t a, b;
        fq_default_poly_t fq_poly;

        // Initalise fq_variables
        fq_default_init(a, ctx);
        fq_default_init(b, ctx);

        // Random polynomial
        fq_default_poly_init(fq_poly, ctx);
        fq_default_poly_randtest(fq_poly, state, n_randint(state, 10), ctx);

        // Evaluate the polynomial twice and ensure results match
        fq_default_randtest(a, state, ctx);
        fq_default_poly_evaluate_fq_default(b, fq_poly, a, ctx);
        fq_default_poly_evaluate_fq_default(a, fq_poly, a, ctx);

        result = (fq_default_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fq_default_poly_print(fq_poly, ctx), flint_printf("\n\n");
            fq_default_print(a, ctx), flint_printf("\n\n");
            fq_default_print(b, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        // Clear everything
        fq_default_poly_clear(fq_poly, ctx);
        fq_default_clear(a, ctx);
        fq_default_clear(b, ctx);
    }
}

TEST_FUNCTION_START(fq_default_poly_evaluate, state)
{
    fmpz_t p;
    fq_default_ctx_t ctx;

    fmpz_init(p);

    // Test FQ_ZECH type with GF(5^5)
    fmpz_set_ui(p, 5);
    fq_default_ctx_init_type(ctx, p, 5, "x", 1);
    test_eval_with_ctx(ctx, state);
    fq_default_ctx_clear(ctx);

    // Test FQ_NMOD type with GF(163^3)
    fmpz_set_ui(p, 163);
    fq_default_ctx_init_type(ctx, p, 3, "x", 2);
    test_eval_with_ctx(ctx, state);
    fq_default_ctx_clear(ctx);

    // Test FQ type with GF((2^127 - 1)^2)
    fmpz_set_str(p, "170141183460469231731687303715884105727", 10);
    fq_default_ctx_init_type(ctx, p, 2, "x", 3);
    test_eval_with_ctx(ctx, state);
    fq_default_ctx_clear(ctx);

    // Test NMOD type with GF(65537)
    fmpz_set_ui(p, 65537);
    fq_default_ctx_init_type(ctx, p, 1, "x", 4);
    test_eval_with_ctx(ctx, state);
    fq_default_ctx_clear(ctx);

    // Test FMPZ_MOD type with GF(2^127 - 1)
    fmpz_set_str(p, "170141183460469231731687303715884105727", 10);
    fq_default_ctx_init_type(ctx, p, 1, "x", 5);
    test_eval_with_ctx(ctx, state);
    fq_default_ctx_clear(ctx);

    fmpz_clear(p);

    TEST_FUNCTION_END(state);
}
