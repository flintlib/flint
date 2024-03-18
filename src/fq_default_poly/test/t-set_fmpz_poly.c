/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fq_default_poly.h"

TEST_FUNCTION_START(fq_default_poly_set_fmpz_poly, state)
{
    int i, result;

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_poly_t fq_poly1, fq_poly2;
        fmpz_poly_t poly;
        fmpz_t p;

        fmpz_init(p);

        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 3, "x");

        fq_default_poly_init(fq_poly1, ctx);
        fq_default_poly_init(fq_poly2, ctx);

        fmpz_poly_init(poly);

        fq_default_poly_one(fq_poly1, ctx);

        fmpz_poly_one(poly);
        fq_default_poly_set_fmpz_poly(fq_poly2, poly, ctx);

        result = fq_default_poly_equal(fq_poly1, fq_poly2, ctx);
        if (!result)
        {
            flint_printf("Polynomials not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(poly);

        fq_default_poly_clear(fq_poly1, ctx);
        fq_default_poly_clear(fq_poly2, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_poly_t fq_poly1, fq_poly2;
        fmpz_poly_t poly;
        fmpz_t p;

        fmpz_init(p);

        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 16, "x");

        fq_default_poly_init(fq_poly1, ctx);
        fq_default_poly_init(fq_poly2, ctx);

        fmpz_poly_init(poly);

        fq_default_poly_one(fq_poly1, ctx);

        fmpz_poly_one(poly);
        fq_default_poly_set_fmpz_poly(fq_poly2, poly, ctx);

        result = fq_default_poly_equal(fq_poly1, fq_poly2, ctx);
        if (!result)
        {
            flint_printf("Polynomials not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(poly);

        fq_default_poly_clear(fq_poly1, ctx);
        fq_default_poly_clear(fq_poly2, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_poly_t fq_poly1, fq_poly2;
        fmpz_poly_t poly;
        fmpz_t p;

        fmpz_init(p);

        fmpz_set_str(p, "73786976294838206473", 10);

        fq_default_ctx_init(ctx, p, 1, "x");

        fq_default_poly_init(fq_poly1, ctx);
        fq_default_poly_init(fq_poly2, ctx);

        fmpz_poly_init(poly);

        fq_default_poly_one(fq_poly1, ctx);

        fmpz_poly_one(poly);
        fq_default_poly_set_fmpz_poly(fq_poly2, poly, ctx);

        result = fq_default_poly_equal(fq_poly1, fq_poly2, ctx);
        if (!result)
        {
            flint_printf("Polynomials not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(poly);

        fq_default_poly_clear(fq_poly1, ctx);
        fq_default_poly_clear(fq_poly2, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
