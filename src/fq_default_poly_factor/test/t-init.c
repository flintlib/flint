/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default_poly_factor.h"

TEST_FUNCTION_START(fq_default_poly_factor_init, state)
{
    int i;

    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_poly_factor_t fq_poly_fac;
        fq_default_poly_t poly;
        fq_default_t c;
        fmpz_t p;

        fmpz_init(p);

        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 3, "x");

        fq_default_poly_factor_init(fq_poly_fac, ctx);

        fq_default_poly_init(poly, ctx);

        fq_default_init(c, ctx);

        while (fq_default_poly_is_zero(poly, ctx))
            fq_default_poly_randtest(poly, state, n_randint(state, 10), ctx);

        fq_default_poly_factor(fq_poly_fac, c, poly, ctx);

        fq_default_clear(c, ctx);

        fq_default_poly_factor_clear(fq_poly_fac, ctx);

        fq_default_ctx_clear(ctx);

        fq_default_ctx_init(ctx, p, 16, "x");

        fq_default_poly_factor_init(fq_poly_fac, ctx);

        fq_default_poly_init(poly, ctx);

        fq_default_init(c, ctx);

        while (fq_default_poly_is_zero(poly, ctx))
            fq_default_poly_randtest(poly, state, n_randint(state, 10), ctx);

        fq_default_poly_factor(fq_poly_fac, c, poly, ctx);

        fq_default_clear(c, ctx);

        fq_default_poly_factor_clear(fq_poly_fac, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_set_str(p, "73786976294838206473", 10);

        fq_default_ctx_init(ctx, p, 1, "x");

        fq_default_poly_factor_init(fq_poly_fac, ctx);

        fq_default_poly_init(poly, ctx);

        fq_default_init(c, ctx);

        while (fq_default_poly_is_zero(poly, ctx))
            fq_default_poly_randtest(poly, state, n_randint(state, 10), ctx);

        fq_default_poly_factor(fq_poly_fac, c, poly, ctx);

        fq_default_clear(c, ctx);

        fq_default_poly_factor_clear(fq_poly_fac, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);

    }

    TEST_FUNCTION_END(state);
}
