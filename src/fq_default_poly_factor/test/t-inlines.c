/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default_poly_factor.h"

TEST_FUNCTION_START(fq_default_poly_factor_inlines, state)
{
    slong i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong deg;
        fq_default_ctx_t ctx;
        fq_default_t t;
        fq_default_poly_t a, b, c, d, e;
        fq_default_poly_factor_t f;

        fmpz_init(p);
        fmpz_randprime(p, state, n_randint(state, 2) ? 40 : 80, 1);
        deg = n_randint(state, 3) + 1;
        fq_default_ctx_init(ctx, p, deg, "a");

        fq_default_init(t, ctx);
        fq_default_poly_init(a, ctx);
        fq_default_poly_init(b, ctx);
        fq_default_poly_init(c, ctx);
        fq_default_poly_init(d, ctx);
        fq_default_poly_init(e, ctx);
        fq_default_poly_factor_init(f, ctx);

        fq_default_poly_one(a, ctx);
        for (j = n_randint(state, 5); j > 0; j--)
        {
            fq_default_poly_randtest_not_zero(b, state, 5, ctx);
            fq_default_poly_mul(a, a, b, ctx);
        }

        fq_default_poly_factor(f, t, a, ctx);

        fq_default_poly_set_fq_default(b, t, ctx);
        for (j = fq_default_poly_factor_length(f, ctx) - 1; j >= 0; j--)
        {
            fq_default_poly_factor_get_poly(c, f, j, ctx);
            FLINT_TEST(fq_default_poly_is_irreducible(c, ctx));
            fq_default_poly_pow(d, c, fq_default_poly_factor_exp(f, j, ctx), ctx);
            fq_default_poly_mul(b, b, d, ctx);
        }

        FLINT_TEST(fq_default_poly_equal(a, b, ctx));

        fq_default_clear(t, ctx);
        fq_default_poly_clear(a, ctx);
        fq_default_poly_clear(b, ctx);
        fq_default_poly_clear(c, ctx);
        fq_default_poly_clear(d, ctx);
        fq_default_poly_clear(e, ctx);
        fq_default_poly_factor_clear(f, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
