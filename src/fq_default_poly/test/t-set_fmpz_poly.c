/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_default_poly.h"

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("set_fmpz_poly....");
    fflush(stdout);

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

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
